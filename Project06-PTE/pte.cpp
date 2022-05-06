#include <fmt/core.h>
#include <yaml-cpp/yaml.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <boost/math/quadrature/gauss.hpp>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <vector>

#include "NumpySaver.hpp"
#include "Timer.hpp"
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
                        int l, int z,
                        std::function<double(double)>& additional_pot) {
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
      double h_val = -values_dd[j - (idx - order)] / 2. +
                     (l * (l + 1.) / (t * t * 2.) - z / t + additional_pot(t)) *
                         values[j - (idx - order)];
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
MatrixXd matrix_H(std::vector<double>& knots, Index order, int l, int z,
                  std::function<double(double)>& additional_pot) {
  Index n = knots.size();
  MatrixXd mat(n - 3 + order, n - 3 + order);

  for (int row = 0; row < mat.rows(); row++)
    for (int col = 0; col < mat.cols(); col++)
      mat(row, col) = matrix_H_element<num_points>(
          row - order + 1, col - order + 1, knots, order, l, z, additional_pot);

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

SparseMatrix<double> matrix_poisson(const std::vector<double>& knots,
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

double rho_hydrogen_ground_state(double r, double e, double r_bohr = 1) {
  return e / M_PI / (r_bohr * r_bohr * r_bohr) * std::exp(-2 * r / r_bohr);
}

VectorXd rhs_poisson(const std::vector<double>& knots,
                     std::function<double(double)>& rho) {
  auto n = knots.size();
  VectorXd rhs(n + 1);
  for (std::size_t i = 0; i < n; i++) {
    rhs[i] = -knots[i] * 4 * M_PI * rho(knots[i]);
  }

  rhs[n] = 0;

  return rhs;
}

double evaluate_spline_poisson(const double& x,
                               const std::vector<double>& knots,
                               const VectorXd& weights, std::size_t order = 3) {
  assert(int(knots.size() + order) - 1 == weights.size());

  auto [values, index] = bsplines::bsplines(x, knots, order);

  double sum = 0;
  for (std::size_t i = 0; i <= order; i++)
    sum += values[i] * weights[index + i];

  return sum;
}

double test_norm_r2(std::function<double(double)> f, double xmin, double xmax) {
  using boost::math::quadrature::gauss;

  std::function<double(double)> integrant = [&](double x) {
    return f(x) * x * x;
  };

  double norm =
      4 * M_PI * gauss<double, 10000>::integrate(integrant, xmin, xmax);

  std::cout << "Norm is " << norm << endl;

  return norm;
}

double test_norm_sqrt(std::function<double(double)> f, double xmin,
                      double xmax) {
  using boost::math::quadrature::gauss;
  std::function<double(double)> integrant = [&](double x) {
    double c = f(x);
    return c * c;
  };

  double norm = gauss<double, 10000>::integrate(integrant, xmin, xmax);

  std::cout << "WVF Norm is " << norm << endl;

  return norm;
}

VectorXd evaluate_spline_poisson(const VectorXd& x,
                                 const std::vector<double>& knots,
                                 const VectorXd& weights,
                                 std::size_t order = 3) {
  auto n = x.size();
  VectorXd ret(n);
  for (Index i = 0; i < n; i++) {
    ret[i] = evaluate_spline(x[i], knots, weights, order);
  }
  return ret;
}

template <int order>
std::vector<double> normalization_constants_radial_wavefunctions(
    std::vector<double>& knots, MatrixXd& eigenvectors) {
  using boost::math::quadrature::gauss;

  std::vector<double> returner;
  returner.reserve(eigenvectors.cols());

  int col;
  std::function<double(double)> f = [&](double r) {
    double v = evaluate_spline(r, knots, eigenvectors.col(col), order);
    return v * v;
  };

  for (col = 0; col < eigenvectors.cols(); col++) {
    double norm = 0;
    for (std::size_t k = 0; k < knots.size() - 1; k++)
      norm += gauss<double, 2 * order>::integrate(f, knots[k], knots[k + 1]);
    returner.push_back(std::sqrt(norm));
  }

  return returner;
}

// int, int, int = N, n, l
template <int order>
std::tuple<double, std::size_t> solve(
    std::vector<double>& knots_atom, std::vector<double>& knots_poisson,
    std::vector<std::tuple<int, int, int>>& orbitals, int z,
    std::string& basefilename, double mixture = 0.4, double tol = 1e-6,
    std::size_t maxiter = 100) {
  using namespace std;
  using boost::math::quadrature::gauss;

  constexpr int gauss_points = double(order) * 1.5;

  size_t n_orbitals = orbitals.size();

  auto matrix_poisson_ = matrix_poisson(knots_poisson);
  // first do the (sparse) LU decomposition
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      solver_poisson;
  solver_poisson.analyzePattern(matrix_poisson_);
  solver_poisson.factorize(matrix_poisson_);

  auto B = matrix_B<gauss_points>(knots_atom, order);

  GeneralizedEigenSolver<MatrixXd> ges;

  VectorXd old_orbital_energies = VectorXd::Zero(n_orbitals);
  VectorXd orbital_energies = VectorXd::Ones(n_orbitals);

  VectorXd weights_poisson = VectorXd::Zero(matrix_poisson_.cols() + 1);
  VectorXd weights_poisson_old = VectorXd::Zero(matrix_poisson_.cols() + 1);

  std::map<int, std::pair<MatrixXd, VectorXd>> l_map, l_map_old;

  std::function rho = [&](double r) -> double {
    if (r == 0) return 0;
    double sum = 0;
    for (auto& [N, n, l] : orbitals) {
      if (!l_map.contains(l)) continue;
      double pnl_r =
          evaluate_spline(r, knots_atom, l_map[l].first.col(n - l - 1), order) /
          r;
      sum += N * pnl_r * pnl_r;
    }
    return sum / 4. / M_PI;
  };

  std::function rho_old = [&](double r) -> double {
    if (r == 0) return 0;
    double sum = 0;
    for (auto& [N, n, l] : orbitals) {
      if (!l_map_old.contains(l)) continue;
      double pnl_r = evaluate_spline(r, knots_atom,
                                     l_map_old[l].first.col(n - l - 1), order) /
                     r;
      sum += N * pnl_r * pnl_r;
    }
    return sum / 4. / M_PI;
  };

  std::function<double(double)> additional = [&](double r) -> double {
    if (r == 0.0) return 0;
    double old_exchange = -3. * pow(3. * rho_old(r) / 8. / M_PI, 1. / 3.);
    double new_exchange = -3. * pow(3. * rho(r) / 8. / M_PI, 1. / 3.);
    double old_direct =
        evaluate_spline_poisson(r, knots_poisson, weights_poisson_old, order) /
        r;
    double new_direct =
        evaluate_spline_poisson(r, knots_poisson, weights_poisson, order) / r;

    return (1. - mixture) * (new_direct + new_exchange) +
           mixture * (old_direct + old_exchange);
  };

  double energy = numeric_limits<double>::max();

  std::size_t iter = 0;

  while ((old_orbital_energies - orbital_energies).norm() > tol) {
    iter++;
    old_orbital_energies = orbital_energies;

    std::map<int, std::pair<MatrixXd, VectorXd>> l_map_new;

    for (size_t i = 0; i < n_orbitals; i++) {
      auto& [N, n, l] = orbitals[i];

      if (!l_map_new.contains(l)) {
        auto H = matrix_H<gauss_points>(knots_atom, order, l, z, additional);

        ges.compute(H, B);
        VectorXcd eigenvalues = ges.eigenvalues();

#ifndef NDEBUG
        if (eigenvalues.imag().cwiseAbs().sum() > 10e-3) {
          cerr << "Error: Complex eigenvalues!" << endl;
        }
#endif

        VectorXd real_eigen = eigenvalues.real();

        // sort the eigenvalues (and create an index list)
        // https://stackoverflow.com/questions/39693909/sort-eigen-matrix-column-values-by-ascending-order-of-column-1-values
        std::vector<int> indices(real_eigen.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&](int A, int B) -> bool {
          return real_eigen[A] < real_eigen[B];
        });

        std::sort(real_eigen.data(), real_eigen.data() + real_eigen.size());

        MatrixXd unsorted_eigenvectors = ges.eigenvectors().real();
        l_map_new[l] = std::pair(MatrixXd(unsorted_eigenvectors.rows(),
                                          unsorted_eigenvectors.cols()),
                                 real_eigen);

        auto normalization_constants =
            normalization_constants_radial_wavefunctions<order>(
                knots_atom, unsorted_eigenvectors);

        for (Index row = 0; row < unsorted_eigenvectors.rows(); row++)
          for (Index col = 0; col < unsorted_eigenvectors.cols(); col++)
            l_map_new[l].first(row, col) =
                unsorted_eigenvectors(row, indices[col]) /
                normalization_constants[indices[col]];

        auto test = normalization_constants_radial_wavefunctions<order>(
            knots_atom, l_map_new[l].first);
      }
      orbital_energies[i] = l_map_new[l].second[n - l - 1];
    }

    l_map_old = l_map;
    l_map = l_map_new;

    auto rhs = rhs_poisson(knots_poisson, rho);

    weights_poisson_old = weights_poisson;
    weights_poisson << 0, solver_poisson.solve(rhs);

    if (iter >= maxiter) {
      cerr << "Reached maximum number of iterations!" << endl;
      break;
    }
  }
  energy = 0;

  function<double(double, int, int)> ee_interaction = [&](double x, int n,
                                                          int l) {
    double pnl_r =
        evaluate_spline(x, knots_atom, l_map[l].first.col(n - l - 1), order);
    return additional(x) * pnl_r * pnl_r;
  };

  for (size_t i = 0; i < n_orbitals; i++) {
    auto& [N, n, l] = orbitals[i];
    function<double(double)> integrant =
        bind(ee_interaction, std::placeholders::_1, n, l);

    double integral = 0;
    for (size_t k = 0; k < knots_poisson.size() - 1; k++) {
      integral += gauss<double, 2 * order>::integrate(
          integrant, knots_poisson[k], knots_poisson[k + 1]);
    }
    energy += N * (orbital_energies[i] - .5 * integral);
  }
  cout << "Orbital energies=" << orbital_energies.transpose() << endl;
  cout << "Total energy = " << energy << endl;

  // exporting/plotting
  VectorXd xx = VectorXd::LinSpaced(1000, 0, 10);
  VectorXd rho_x = VectorXd(xx.size());
  for (Index i = 0; i < xx.size(); i++) rho_x[i] = rho(xx[i]);
  NumpySaver(fmt::format("build/output/{}_rho.npy", basefilename))
      << xx << rho_x;

  return {energy, iter};
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

std::vector<std::tuple<int, int, int>> to_tuples(
    std::vector<std::vector<int>>& list) {
  // how i wish map was a thing here...

  std::vector<std::tuple<int, int, int>> orbitals;
  orbitals.reserve(list.size());
  for (auto& orbit : list)
    orbitals.push_back(std::make_tuple(orbit[0], orbit[1], orbit[2]));
  return orbitals;
}

int main() {
  int n_knots = 50;
  int n_knots_p = 1000;

  std::vector<double> knots_atom(n_knots);
  std::vector<double> knots_poisson(n_knots_p);

  for (int i = 0; i < n_knots; i++)
    knots_atom[i] = (1000. - 0.) *
                        (std::exp(double(i) / double(n_knots - 1) * 10.) - 1.) /
                        (std::exp(10.) - 1.) +
                    0.;
  for (int i = 0; i < n_knots_p; i++)
    knots_poisson[i] =
        (1000. - 0.) *
            (std::exp(double(i) / double(n_knots_p - 1) * 10.) - 1.) /
            (std::exp(10.) - 1.) +
        0.;

  Map<VectorXd> knots_map(knots_atom.data(), knots_atom.size());
  Map<VectorXd> knots_map_p(knots_poisson.data(), knots_poisson.size());
  NumpySaver("build/output/knots_poisson.npy") << knots_map_p;
  NumpySaver("build/output/knots_atom.npy") << knots_map;

  // Check the influence of the mixing factor
  std::vector<std::tuple<int, int, int>> orbitals_Ne = {
      {2, 1, 0}, {2, 2, 0}, {6, 2, 1}};  // N, n, l
  std::vector<std::tuple<int, int, int>> orbitals_Nep = {
      {2, 1, 0}, {2, 2, 0}, {5, 2, 1}};  // N, n, l

  VectorXd etas = VectorXd::LinSpaced(50, 0, 1);
  VectorXd iters(etas.size());
  for (Index i = 0; i < etas.size(); i++) {
    std::string name = fmt::format("eta_{}", etas[i]);
    auto [energy, iter] =
        solve<3>(knots_atom, knots_poisson, orbitals_Ne, 10, name, etas[i]);
    iters[i] = iter;
  }
  NumpySaver("build/output/etas.npy") << etas << iters;

  // iterate through ALL the elements
  YAML::Node elements = YAML::LoadFile("build/output/pse.yaml");

  std::size_t num_elements = elements.size();

  VectorXd zs(num_elements);
  VectorXd es(num_elements);
  VectorXd iterss(num_elements);
  VectorXd es_ion(num_elements);
  VectorXd iterss_ion(num_elements);

  for (std::size_t i = 0; i < elements.size(); i++) {
    YAML::Node element = elements[i];

    auto z = element["Z"].as<int>();
    auto name = element["name"].as<std::string>();
    auto symbol = element["symbol"].as<std::string>();
    auto orbitals_ = element["orbitals"].as<std::vector<std::vector<int>>>();
    auto orbitals_ionized_ =
        element["orbitals_ionized"].as<std::vector<std::vector<int>>>();

    auto orbitals = to_tuples(orbitals_);
    auto orbitals_ionized = to_tuples(orbitals_ionized_);

    cout << "Element " << name << "(" << i + 1 << "/" << elements.size() << ")"
         << endl;
    zs[i] = z;

    Timer t;
    auto [energy, iter] =
        solve<3>(knots_atom, knots_poisson, orbitals, z, symbol, .6);
    cout << "E=" << energy << " took " << iter << " iterations ("
         << t.elapsed_s() << "s)" << endl;

    es[i] = energy;
    iterss[i] = iter;

    std::string symbol_ionized = symbol + "+";

    t.reset();
    auto [energy_ion, iter_ion] = solve<3>(
        knots_atom, knots_poisson, orbitals_ionized, z, symbol_ionized, .4);
    cout << "E_{ion}=" << energy_ion << " took " << iter_ion << " iterations ("
         << t.elapsed_s() << "s)" << endl;

    es_ion[i] = energy_ion;
    iterss_ion[i] = iter_ion;
    cout << endl;
  }

  NumpySaver("build/output/elements.npy")
      << zs << es << iterss << es_ion << iterss_ion;
}
