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
  double r =
      std::sqrt((x[0] - x[2]) * (x[0] - x[2]) + (x[1] - x[3]) * (x[1] - x[3]));
  return std::exp(-x[0] * x[0] - x[1] * x[1] - x[2] * x[2] - x[3] * x[3] +
                  2. * lambda * r / (1. + alpha * r));
}

constexpr double pow2(double x) { return x * x; };
constexpr double pow3(double x) { return x * x * x; };

constexpr double derivative_factor(double r, double x, double x_other,
                                   double alpha, double lambda) {
  using std::pow;
  return -1. + lambda / pow2(1 + alpha * r) / r -
         lambda * pow2(x - x_other) / pow3((1. + alpha * r) * r) *
             (1. + 3. * alpha * r) +
         pow2(-x + lambda * (x - x_other) / pow2(1. + alpha * r) / r);
}

double energy_2d(const Array4d& x, double alpha, double lambda) {
  double r =
      std::sqrt((x[0] - x[2]) * (x[0] - x[2]) + (x[1] - x[3]) * (x[1] - x[3]));

  double dx1 = derivative_factor(r, x[0], x[2], alpha, lambda);
  double dx2 = derivative_factor(r, x[2], x[1], alpha, lambda);
  double dy1 = derivative_factor(r, x[1], x[3], alpha, lambda);
  double dy2 = derivative_factor(r, x[3], x[1], alpha, lambda);

  return (-dx1 + pow2(x[0]) - dx2 + pow2(x[2]) - dy1 + pow2(x[1]) - dy2 +
          pow2(x[3])) /
             2. +
         lambda / r;
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

  NumpySaver("build/output/test_2D_distribution.npy") << test_1;
}

void test_1D_energy() {
  Array1d argstart;
  argstart << 0;

  std::normal_distribution<> step_gauss(0, 1);

  using namespace std::placeholders;
  auto pdf_normal = std::bind(pdf_1d, _1, 1. / 4.);

  MetropolisAlgorithm<1, std::normal_distribution<> > malg_normal(
      pdf_normal, argstart, step_gauss);
}
int main(int argc, char const* argv[]) {
  test_1D_metropolis();
  test_2D_metropolis();
  return 0;
}
