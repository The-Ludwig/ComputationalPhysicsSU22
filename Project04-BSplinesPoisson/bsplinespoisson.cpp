
#include <fmt/core.h>

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cassert>
#include <iostream>
#include <vector>

#include "NumpySaver.hpp"
#include "Timer.hpp"
#include "YAMLUtils.hpp"
#include "bsplines.hpp"

using namespace Eigen;
using std::cout, std::endl, std::cerr;

void test_bsplines(double x = .1, Index order = 3, Index nknots = 11) {
  std::vector<double> std_knots(nknots);
  Map<ArrayXd> knots(std_knots.data(), std_knots.size());
  knots = ArrayXd::LinSpaced(nknots, 0, 1);

  auto [values, index] = bsplines::bsplines(x, std_knots, order);

  Map<VectorXd> values_eigen(values.data(), values.size());
  auto sum = values_eigen.sum();
  if (std::abs(sum - 1) > 1e-15) {
    cerr << "The knots are: " << knots.transpose() << endl;
    cerr << "The spline value at " << x << " is: " << values_eigen.transpose()
         << endl;
    cerr << "Summing all spline values (should be 1): " << sum << endl;
  }
}

SparseMatrix<double> build_matrix(const std::vector<double>& knots,
                                  std::size_t order = 3) {
  auto n = knots.size();

  SparseMatrix<double> matrix(n + 1, n + order - 2);
  matrix.reserve(order * (n + 1) - 2);

  auto [deriv_0, _index] = bsplines::ndx_bsplines(knots[0], knots, order, 2);
  for (int col = 0; col < int(order) - 1; col++)
    matrix.insert(0, col) = deriv_0[col + 1];

  for (int row = 1; row < int(n); row++) {
    auto [deriv, index] = bsplines::ndx_bsplines(knots[row], knots, order, 2);
    for (int col = row - 1; col < row - 1 + int(order); col++)
      matrix.insert(row, col) = deriv[col - row + 2];
  }

  auto [deriv_last, __index] =
      bsplines::ndx_bsplines(knots[n - 1], knots, order, 1);
  matrix.insert(n, n + order - 3) = deriv_last[order];
  matrix.insert(n, n + order - 4) = deriv_last[order - 1];

  return matrix;
}

double rho_uniform_sphere(double r, double R, double q) {
  if (r <= R) return q / (4. * M_PI / 3. * R * R * R);

  return 0;
}

double rho_uniform_shell(double r, double R1, double R2, double q) {
  if (r < R1) return 0;
  if (r <= R2) return q / (4. * M_PI / 3. * (R2 * R2 * R2 - R1 * R1 * R1));

  return 0;
}

double rho_hydrogen_ground_state(double r, double e, double r_bohr = 1) {
  return e / M_PI / (r_bohr * r_bohr * r_bohr) * std::exp(-2 * r / r_bohr);
}

VectorXd build_rhs(const std::vector<double>& knots,
                   std::function<double(double)>& rho) {
  auto n = knots.size();
  VectorXd rhs(n + 1);
  for (std::size_t i = 0; i < n; i++) {
    rhs[i] = -knots[i] * 4 * M_PI * rho(knots[i]);
  }

  rhs[n] = 0;

  return rhs;
}

double evaluate_spline(const double& x, const std::vector<double>& knots,
                       const VectorXd& weights, std::size_t order = 3) {
  auto n = knots.size();
  assert(int(n + order) - 1 == weights.size());

  auto [values, index] = bsplines::bsplines(x, knots, order);

  double sum = 0;
  for (std::size_t i = 0; i <= order; i++)
    sum += values[i] * weights[index + i];

  return sum;
}

VectorXd evaluate_spline(const VectorXd& x, const std::vector<double>& knots,
                         const VectorXd& weights, std::size_t order = 3) {
  auto n = x.size();
  VectorXd ret(n);
  for (Index i = 0; i < n; i++) {
    ret[i] = evaluate_spline(x[i], knots, weights, order);
  }
  return ret;
}

void solve_problem(YAML::Node node) {
  using namespace std::placeholders;

  auto name = node["name"].as<std::string>();
  auto R = node["R"].as<double>();
  auto R1 = node["R1"].as<double>();
  auto R2 = node["R2"].as<double>();
  auto q = node["q"].as<double>();
  auto e = node["e"].as<double>();
  auto r_bohr = node["r_bohr"].as<double>();

  YAML::Node knot_node = node["knots"];
  auto knots = get_yaml_values<double>(knot_node);

  auto matrix = build_matrix(knots);

  std::vector<std::tuple<const char*, std::function<double(double)>>> rhos = {
      {"uniform_sphere", std::bind(rho_uniform_sphere, _1, R, q)},
      {"uniform_shell", std::bind(rho_uniform_shell, _1, R1, R2, q)},
      {"hydrogen", std::bind(rho_hydrogen_ground_state, _1, e, r_bohr)},
  };

  if (knots.size() < 30) {
    cout << "Knots = ";
    for (std::size_t i = 0; i < knots.size(); i++) cout << knots[i] << ", ";
    cout << "\n";
  }

  // first do the (sparse) LU decomposition
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      solver;
  solver.analyzePattern(matrix);
  solver.factorize(matrix);

  VectorXd x = VectorXd::LinSpaced(1000, knots[0], knots[knots.size() - 1]);

  for (std::size_t i = 0; i < rhos.size(); i++) {
    auto rho_name = std::get<0>(rhos[i]);
    auto rho = std::get<1>(rhos[i]);

    cout << "\nProblem = " << rho_name << "\n";

    auto rhs = build_rhs(knots, rho);
    auto solution = solver.solve(rhs);
    if (knots.size() < 30) {
      cout << "rhs vector: " << rhs.transpose() << "\n";
      cout << "solution vector: " << solution.transpose() << "\n";
    }

    VectorXd weights(solution.size() + 1);
    weights << 0, solution;

    VectorXd sol = evaluate_spline(x, knots, weights);

    NumpySaver(fmt::format("build/output/solution_{}_{}.npy", name, rho_name))
        << x << sol;
  }
  cout << "\n";
}

void plot_bsplines(std::size_t order = 3) {
  std::vector<double> std_knots(11);
  Map<ArrayXd> knots(std_knots.data(), std_knots.size());
  knots = ArrayXd::LinSpaced(11, 0, 1);
  VectorXd xs = VectorXd::LinSpaced(1000, knots[0], knots[knots.size() - 1]);
  MatrixXd spline(xs.size(), order + 1);
  MatrixXd dspline(xs.size(), order + 1);
  MatrixXd ddspline(xs.size(), order + 1);

  for (Index i = 0; i < xs.size(); i++) {
    auto [x_spline, index_0] = bsplines::bsplines(xs[i], std_knots, order);
    auto [x_dspline, index_1] =
        bsplines::ndx_bsplines(xs[i], std_knots, order, 1);
    auto [x_ddspline, index_2] =
        bsplines::ndx_bsplines(xs[i], std_knots, order, 2);
    for (std::size_t j = 0; j < order + 1; j++) {
      spline(i, j) = x_spline[j];
      dspline(i, j) = x_dspline[j];
      ddspline(i, j) = x_ddspline[j];
    }
  }
  NumpySaver saver("build/output/plot_splines.npy");
  saver << xs;
  for (std::size_t j = 0; j < order + 1; j++) saver << spline.col(j);
  for (std::size_t j = 0; j < order + 1; j++) saver << dspline.col(j);
  for (std::size_t j = 0; j < order + 1; j++) saver << ddspline.col(j);
}

int main() {
  cout << "Testing if all spline values at each point sum to 1, no error means "
          "it worked."
       << endl;
  for (double x = 0; x <= 1; x += .001) {
    test_bsplines(x, 3, 11);
  }

  std::vector<double> std_knots(11);
  Map<ArrayXd> knots(std_knots.data(), std_knots.size());
  knots = ArrayXd::LinSpaced(11, 0, 1);

  cout << "Testing Matrix creation..." << endl;
  auto matrix = build_matrix(std_knots);
  cout << "Matrix size is " << matrix.rows() << "x" << matrix.cols() << endl;
  cout << matrix << endl;

  YAML::Node config = YAML::LoadFile("config.yaml");

  for (std::size_t i = 0; i < config.size(); i++) {
    cout << "\n\n*****************************************\n"
         << fmt::format("Running simulation {}/{}", i + 1, config.size())
         << endl;
    solve_problem(config[i]);
  }

  plot_bsplines();
  return 0;
}
