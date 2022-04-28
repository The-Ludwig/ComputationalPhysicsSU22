#include <fmt/core.h>
#include <yaml-cpp/yaml.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <boost/math/quadrature/gauss.hpp>
#include <iostream>
#include <numeric>
#include <optional>
#include <vector>

#include "NumpySaver.hpp"
#include "bsplines.hpp"

using namespace Eigen;
using std::cout, std::endl, std::cerr;

const char* devider = "##########################################";

template <int num_points>
double matrix_B_element(int i, int j, std::vector<double>& knots, Index order) {
  using boost::math::quadrature::gauss;

  if (std::abs((int)i - (int)j) > order) return 0;

  int imax = std::max(i, j);
  int imin = std::min(i, j);

  Index n = knots.size();

  double sum = 0;

  for (int k = std::max(0, imax - 1); k < std::min(imin + order + 1, n - 1);
       k++) {
    auto f = [&](const double& t) {
      auto [values, idx] = bsplines::bsplines(t, knots, order);
      return values[i + order - idx] * values[j + order - idx];
    };
    sum += gauss<double, num_points>::integrate(f, knots[k], knots[k + 1]);
  }

  return sum;
}

template <int num_points>
double matrix_H_element(int i, int j, std::vector<double>& knots, Index order,
                        int l, int z) {
  using boost::math::quadrature::gauss;
  if (std::abs((int)i - (int)j) > order) return 0;

  Index n = knots.size();
  int imax = std::max(i, j);
  int imin = std::min(i, j);

  double sum = 0;

  for (int k = std::max(0, imax - 1); k < std::min(imin + order + 1, n - 1);
       k++) {
    auto f = [&](const double& t) {
      auto [values, idx] = bsplines::bsplines(t, knots, order);
      auto [values_dd, idx_dd] = bsplines::ndx_bsplines(t, knots, order, 2);
      double h_val =
          -values_dd[j - (idx - order)] / 2. +
          (l * (l + 1.) / (t * t * 2.) - z / t) * values[j - (idx - order)];
      return values[i - (idx - order)] * h_val;
      ;
    };

    sum += gauss<double, num_points>::integrate(f, knots[k], knots[k + 1]);
  }
  return sum;
}

template <int num_points>
MatrixXd matrix_B(std::vector<double>& knots, Index order) {
  Index n = knots.size();
  MatrixXd mat(n - 3 + order, n - 3 + order);

  for (int row = 0; row < mat.rows(); row++)
    for (int col = 0; col < mat.cols(); col++)
      mat(row, col) = matrix_B_element<num_points>(
          row - order + 1, col - order + 1, knots, order);

  return mat;
}

template <int num_points>
MatrixXd matrix_H(std::vector<double>& knots, Index order, int l, int z) {
  Index n = knots.size();
  MatrixXd mat(n - 3 + order, n - 3 + order);

  for (int row = 0; row < mat.rows(); row++)
    for (int col = 0; col < mat.cols(); col++)
      mat(row, col) = matrix_H_element<num_points>(
          row - order + 1, col - order + 1, knots, order, l, z);

  return mat;
}

double evaluate_spline(const double& x, const std::vector<double>& knots,
                       const VectorXd& weights, std::size_t order) {
  auto n = knots.size();
  assert(int(n + order) - 3 == weights.size());

  auto [values, index] = bsplines::bsplines(x, knots, order);

  double sum = 0;
  for (std::size_t i = std::max(int(index) - 1, 0);
       i <= std::min(index - 1 + order, n + order - 4); i++)
    sum += values[i + 1 - index] * weights[i];

  return sum;
}

template <int order>
void solve(std::vector<double>& knots, std::string& basefilename,
           std::optional<double> plot_max, int l, int z) {
  constexpr int gauss_points = order * 1.5;
  int n_knots = knots.size();

  cout << "knots = ";
  for (auto k : knots) cout << k << ", ";
  cout << endl;

  auto B = matrix_B<gauss_points>(knots, order);
  auto H = matrix_H<gauss_points>(knots, order, l, z);

  if (knots.size() < 20) {
    cout << "B=\n" << B << "\n\n";
    cout << "H=\n" << H << endl;
  }

  GeneralizedEigenSolver<MatrixXd> ges;
  ges.compute(H, B);
  VectorXcd eigenvalues = ges.eigenvalues();
  VectorXd real_eigen = eigenvalues.real();

  // sort the eigenvalues (and create an index list)
  std::vector<int> indices(real_eigen.size());
  std::iota(indices.begin(), indices.end(), 0);

  std::sort(indices.begin(), indices.end(), [&](int A, int B) -> bool {
    return real_eigen[A] < real_eigen[B];
  });
  std::sort(real_eigen.data(), real_eigen.data() + real_eigen.size());

  if (eigenvalues.imag().cwiseAbs().sum() > 10e-10) {
    cerr << "Error: Complex eigenvalues!" << endl;
  }

  cout << "Eigenvalues: " << real_eigen.transpose() << endl;

  MatrixXd eigenvectors = ges.eigenvectors().real();

  NumpySaver(fmt::format("build/output/{}_eigenvalues.npy", basefilename))
      << real_eigen;
  NumpySaver(fmt::format("build/output/{}_knots.npy", basefilename)) << knots;

  NumpySaver eigvec_saver(
      fmt::format("build/output/{}_eigenvectors.npy", basefilename));
  for (int col = 0; col < eigenvectors.cols(); col++)
    eigvec_saver << eigenvectors.col(indices[col]);

  if (auto xmax = plot_max) {
    ArrayXd x = ArrayXd::LinSpaced(1000, 0, *xmax);
    MatrixXd points(x.size(), eigenvectors.cols());
    for (int row = 0; row < x.size(); row++)
      for (int col = 0; col < points.cols(); col++)
        points(row, col) = evaluate_spline(
            x[row], knots, eigenvectors.col(indices[col]), order);

    NumpySaver plot_saver(
        fmt::format("build/output/{}_plot_points.npy", basefilename));
    plot_saver << x;
    for (int col = 0; col < points.cols(); col++) plot_saver << points.col(col);
  }
}

template <int order>
void solve_exp(int n_knots, double rmin, double rmax, double exp_fak,
               std::string& basefilename, std::optional<double> plot, int l,
               int z) {
  std::vector<double> knots(n_knots);
  for (int i = 0; i < n_knots; i++)
    knots[i] = (rmax - rmin) *
                   (std::exp(double(i) / double(n_knots - 1) * exp_fak) - 1.) /
                   (std::exp(exp_fak) - 1) +
               rmin;

  solve<order>(knots, basefilename, plot, l, z);
}

template <int order>
void solve_lin(int n_knots, double rmin, double rmax, std::string& basefilename,
               std::optional<double> plot, int l, int z) {
  std::vector<double> knots(n_knots);
  for (int i = 0; i < n_knots; i++)
    knots[i] = double(i) / double(n_knots - 1) * (rmax - rmin) + rmin;

  solve<order>(knots, basefilename, plot, l, z);
}

int main() {
  YAML::Node configs = YAML::LoadFile("config.yaml");

  constexpr int order = 3;

  for (std::size_t i = 0; i < configs.size(); i++) {
    YAML::Node config = configs[i];

    std::optional<double> plot;

    if (auto p_node = config["plot_xmax"]) {
      try {
        auto val = p_node.as<double>();
        if (val > 0) plot.emplace(val);
      } catch (const YAML::TypedBadConversion<double>&) {
        auto boolval = p_node.as<bool>();
        if (boolval) plot.emplace(10);
      }
    }

    auto l = config["l"].as<int>();
    auto z = config["Z"].as<int>();
    auto basefilename = config["name"].as<std::string>();

    fmt::print("\n\n{0}\n# Simulation '{1}' {2}/{3}\n{0}\n", devider,
               basefilename, i + 1, configs.size());

    YAML::Node knot_node = config["knots"];
    auto type = knot_node.Type();
    // single value
    if (type == YAML::NodeType::Sequence) {
      auto knots = knot_node.as<std::vector<double>>();
      solve<order>(knots, basefilename, plot, l, z);
    } else if (type == YAML::NodeType::Map) {
      auto rmin = knot_node["rmin"].as<double>();
      auto rmax = knot_node["rmax"].as<double>();
      auto N = knot_node["N"].as<size_t>();
      auto type = knot_node["type"].as<std::string>();
      if (type == "linear") {
        solve_lin<order>(N, rmin, rmax, basefilename, plot, l, z);
      } else if (type == "exponential") {
        auto expfac = knot_node["expfac"].as<double>();
        solve_exp<order>(N, rmin, rmax, expfac, basefilename, plot, l, z);
      } else {
        throw std::runtime_error(
            "'type' of knots must either be 'linear' or 'exponential'.");
      }
    } else {
      throw std::runtime_error(
          "Range node must either be list, scalar or a map defining "
          "'type', 'N', 'rmin' and 'rmax'.");
    }
  }
}