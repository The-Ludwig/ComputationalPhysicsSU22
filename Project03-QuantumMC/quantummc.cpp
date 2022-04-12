#include <fmt/core.h>
#include <yaml-cpp/yaml.h>

#include <Eigen/Dense>
#include <iostream>
#include <random>

#include "Metropolis.hpp"
#include "NumpySaver.hpp"
#include "Timer.hpp"

using namespace Eigen;
using std::cout, std::endl;

using Array1d = Array<double, 1, 1>;

constexpr double pow2(double x) { return x * x; };
constexpr double pow3(double x) { return x * x * x; };

double pdf_1d(const Array1d& x, double alpha) {
  return std::exp(-2 * alpha * x[0] * x[0]);
}

double energy_1d(const Array1d& x, double alpha) {
  return alpha + x[0] * x[0] * (1. / 2. - 2 * alpha * alpha);
}

/**
 * @brief
 *
 * @param x = {x0, y0, x1, y1}
 * @param alpha
 * @param lambda
 * @return double
 */
double pdf_wavefunction_2d(const Array4d& x, double alpha, double lambda) {
  double r = std::sqrt(pow2(x[0] - x[2]) + pow2(x[1] - x[3]));
  return std::exp(-pow2(x[0]) - pow2(x[1]) - pow2(x[2]) - pow2(x[3]) +
                  2. * lambda * r / (1. + alpha * r));
}

constexpr double derivative_factor(double r, double x, double x_other,
                                   double alpha, double lambda) {
  using std::pow;
  return -1. + lambda / pow2(1 + alpha * r) / r -
         lambda * pow2(x - x_other) / pow3((1. + alpha * r) * r) *
             (1. + 3. * alpha * r) +
         pow2(-x + lambda * (x - x_other) / pow2(1. + alpha * r) / r);
}

double energy_2d(const Array4d& x, double alpha, double lambda) {
  double r = std::sqrt(pow2(x[0] - x[2]) + pow2(x[1] - x[3]));

  double dx1 = derivative_factor(r, x[0], x[2], alpha, lambda);
  double dx2 = derivative_factor(r, x[2], x[0], alpha, lambda);
  double dy1 = derivative_factor(r, x[1], x[3], alpha, lambda);
  double dy2 = derivative_factor(r, x[3], x[1], alpha, lambda);

  return (-dx1 + pow2(x[0]) - dx2 + pow2(x[2]) - dy1 + pow2(x[1]) - dy2 +
          pow2(x[3])) /
             2. +
         lambda / r;
}

std::tuple<double, double> golden_search_simple(std::function<double(double)> f,
                                                double xmin, double xmax,
                                                double tol) {
  constexpr double golden_ratio = 1.618033988749895;
  double x1, x2;
  if (xmin >= xmax)
    throw std::invalid_argument("xmin must be smaller than xmax");

  while (xmax - xmin > tol) {
    x1 = xmax - (xmax - xmin) / golden_ratio;
    x2 = xmin + (xmax - xmin) / golden_ratio;

    if (f(x1) < f(x2))
      xmax = x2;
    else
      xmin = x1;
  }

  return {xmin, xmax};
}

std::tuple<double, double> golden_search(std::function<double(double)> f,
                                         double xmin, double xmax, double tol) {
  constexpr double golden_ratio = 1.618033988749895;
  double x1 = xmax - (xmax - xmin) / golden_ratio;
  double x2 = xmin + (xmax - xmin) / golden_ratio;
  double f1 = f(x1);
  double f2 = f(x2);

  if (xmin >= xmax)
    throw std::invalid_argument("xmin must be smaller than xmax");

  while (xmax - xmin > tol) {
    if (f1 < f2) {
      xmax = x2;
      x2 = x1;
      f2 = f1;
      x1 = xmax - (xmax - xmin) / golden_ratio;
      f1 = f(x1);
    } else {
      xmin = x1;
      x1 = x2;
      f1 = f2;
      x2 = xmin + (xmax - xmin) / golden_ratio;
      f2 = f(x2);
    }
  }

  return {xmin, xmax};
}

std::tuple<double, double> adaptive_golden_search(
    std::function<std::tuple<double, double>(double, Index sample)> f,
    double xmin, double xmax, double tol, Index start_sample = 1000) {
  constexpr double golden_ratio = 1.618033988749895;
  Index sample = start_sample;

  double x1 = xmax - (xmax - xmin) / golden_ratio;
  double x2 = xmin + (xmax - xmin) / golden_ratio;
  auto [f1, f1_err] = f(x1, sample);
  auto [f2, f2_err] = f(x2, sample);

  std::tuple<double, double> buf;

  if (xmin >= xmax)
    throw std::invalid_argument("xmin must be smaller than xmax");

  while (xmax - xmin > tol) {
    if (f1 < f2) {
      xmax = x2;
      x2 = x1;
      f2 = f1;
      x1 = xmax - (xmax - xmin) / golden_ratio;
      std::tie(f1, f1_err) = f(x1, sample);
    } else {
      xmin = x1;
      x1 = x2;
      f1 = f2;
      x2 = xmin + (xmax - xmin) / golden_ratio;
      std::tie(f2, f2_err) = f(x2, sample);
    }
    if (std::abs(f1 - f2) < f1_err + f2_err) {
      sample += (Index)pow2((f1_err + f2_err) / std::abs(f1 - f2));
    }
  }

  return {xmin, xmax};
}

void test_1D_metropolis() {
  Array1d argstart;
  argstart << 0;

  std::uniform_real_distribution<> step(-2, 2);
  std::normal_distribution<> step_gauss(0, 1);

  using namespace std::placeholders;
  auto pdf_normal = std::bind(pdf_1d, _1, 1. / 4.);

  MetropolisAlgorithm<1> malg_uniform(pdf_normal, argstart, step);
  MetropolisAlgorithm<1, std::normal_distribution<> > malg_normal(
      pdf_normal, argstart, step_gauss);

  auto test_1 = malg_uniform.get_sample(100000);
  auto test_2 = malg_normal.get_sample(100000);

  auto [mean, std] =
      Timer::measure_time([&]() { malg_uniform.get_sample(100000); }, 1);
  cout << "Generating 100000 events with uniform step took " << mean << "±"
       << std << "ms\n";

  auto [mean_n, std_n] =
      Timer::measure_time([&]() { malg_normal.get_sample(100000); }, 1);
  cout << "Generating 100000 events with normal step took " << mean_n << "±"
       << std_n << "ms\n";
  NumpySaver("build/output/test_normal_distribution.npy") << test_1 << test_2;
}

void test_2D_metropolis() {
  Array4d argstart;
  argstart << 0, 1, 0, -1;

  std::uniform_real_distribution<> step(-2, 2);
  std::normal_distribution<> step_gauss(0, 1);

  using namespace std::placeholders;
  auto pdf_2d = std::bind(pdf_wavefunction_2d, _1, 1, 1);

  MetropolisAlgorithm<4, std::normal_distribution<> > malg_normal(
      pdf_2d, argstart, step_gauss);

  auto test_1 = malg_normal.get_sample(100000);

  auto [mean, std] =
      Timer::measure_time([&]() { malg_normal.get_sample(100000); }, 1);
  cout << "Generating 100000 2D events with uniform step took " << mean << "±"
       << std << "ms\n";

  NumpySaver("build/output/test_2D_distribution.npy")
      << test_1.col(0) << test_1.col(1) << test_1.col(2) << test_1.col(3);
}

void test_1D_energy(Index n = 300) {
  using namespace std::placeholders;

  Array1d argstart;
  argstart << 0;
  std::normal_distribution<> step_gauss(0, 1);

  ArrayXd alphas = ArrayXd::LinSpaced(n, 0.1, 1);
  ArrayXd means(n);
  ArrayXd stds(n);

  for (Index i = 0; i < n; i++) {
    auto pdf_normal = std::bind(pdf_1d, _1, alphas[i]);
    auto energy = std::bind(energy_1d, _1, alphas[i]);
    MetropolisAlgorithm<1, std::normal_distribution<> > malg_normal(
        pdf_normal, argstart, step_gauss);
    auto [mean, std] = malg_normal.average(energy, 100000);
    means[i] = mean;
    stds[i] = std;
  }

  NumpySaver("build/output/test_1d_energy.npy") << alphas << means << stds;
}

void test_2D_energy(Index n = 300, double lambda = 1) {
  using namespace std::placeholders;

  Array4d argstart;
  argstart << 0, 1, 0, -1;
  std::normal_distribution<> step_gauss(0, 1);

  ArrayXd alphas = ArrayXd::LinSpaced(n, 0.1, 1);
  ArrayXd means(n);
  ArrayXd stds(n);

  for (Index i = 0; i < n; i++) {
    auto pdf = std::bind(pdf_wavefunction_2d, _1, alphas[i], lambda);
    auto energy = std::bind(energy_2d, _1, alphas[i], lambda);
    MetropolisAlgorithm<4, std::normal_distribution<> > malg(pdf, argstart,
                                                             step_gauss);
    auto [mean, std] = malg.average(energy, 100000);
    means[i] = mean;
    stds[i] = std;
  }

  NumpySaver(fmt::format("build/output/test_2d_energy_{}.npy", lambda))
      << alphas << means << stds;
}

void find_2d_alpha(double lambda = 1, Index sample = 1000000) {
  using namespace std::placeholders;

  Array4d argstart;
  argstart << 0, 1, 0, -1;
  std::normal_distribution<> step_gauss(0, 1);

  std::function<double(double)> f = [&](double alpha) {
    auto pdf = std::bind(pdf_wavefunction_2d, _1, alpha, lambda);
    auto energy = std::bind(energy_2d, _1, alpha, lambda);
    MetropolisAlgorithm<4, std::normal_distribution<> > malg(pdf, argstart,
                                                             step_gauss);
    auto [mean, std] = malg.average(energy, sample);
    return mean;
  };

  auto [xmin, xmax] = golden_search(f, .1, 1, 1e-6);
  cout << "Optimal α for 2D case (λ=" << lambda << ") is between " << xmin
       << " and " << xmax << " (" << f(xmin) << " < E <" << f(xmax) << ")"
       << endl;
}

void find_2d_alpha_adaptive(double lambda = 1) {
  using namespace std::placeholders;

  Array4d argstart;
  argstart << 0, 1, 0, -1;
  std::normal_distribution<> step_gauss(0, 1);

  std::function<std::tuple<double, double>(double, Index)> f =
      [&](double alpha, double sample) {
        auto pdf = std::bind(pdf_wavefunction_2d, _1, alpha, lambda);
        auto energy = std::bind(energy_2d, _1, alpha, lambda);
        MetropolisAlgorithm<4, std::normal_distribution<> > malg(pdf, argstart,
                                                                 step_gauss);
        return malg.average(energy, sample);
      };

  auto [xmin, xmax] = adaptive_golden_search(f, .1, 1, 1e-6);
  cout << "Optimal (adaptive) α for 2D case (λ=" << lambda << ") is between "
       << xmin << " and " << xmax << endl;
}

int main(int argc, char const* argv[]) {
  test_1D_metropolis();
  test_2D_metropolis();
  test_1D_energy();
  test_2D_energy(300, 0);
  test_2D_energy(300, 1);
  test_2D_energy(300, 2);
  test_2D_energy(300, 8);
  find_2d_alpha(0);
  find_2d_alpha(1);
  find_2d_alpha(2);
  find_2d_alpha(8);
  return 0;
}
