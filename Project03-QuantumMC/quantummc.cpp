#include <yaml-cpp/yaml.h>

#include <Eigen/Dense>
#include <iostream>
#include <random>

#include "NumpySaver.hpp"
#include "Timer.hpp"

using namespace Eigen;

template <int NumArguments>
Array<double, Dynamic, NumArguments> metropolis(
    std::function<double(const Array<double, NumArguments, 1>&)> pdf,
    unsigned int N, unsigned int starting_phase = 10,
    Array<double, NumArguments, 1> argstart =
        Array<double, NumArguments, 1>::Zero()) {
  typedef Eigen::Array<double, NumArguments, 1> Argument;
  Eigen::ArrayXd returner(N);

  std::random_device
      rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());  // Standard mersenne_twister_engine seeded with rd()

  std::uniform_real_distribution<> uniform(0., 1.0);
  std::uniform_real_distribution<> step(-1., 1.);
  // Here it gets confusing: We are only using one argument
  // We will change the entries of the argument seperately
  // so we don't copy every value unnecesary
  Argument arg = argstart;
  // saving the old p, minimizes the number of function calls
  double p_old = pdf(arg);

  // starting phase (throwing values away)
  for (std::size_t i = 0; i < starting_phase; i++) {
    for (std::size_t j = 0; j < NumArguments; j++) {
      double x_old = arg[j];
      double x_new = x_old + step(gen);

      arg[j] = x_new;
      double p_new = pdf(arg);
      double p = p_new / p_old;

      if (p < 1) {
        double s = uniform(gen);

        // reject new values
        if (s > p) {
          x_new = x_old;
          p_new = p_old;
        }
      }

      arg[j] = x_new;
      p_old = p_new;
    }
  }

  // now saving the values
  for (Eigen::Index i = 0; i < N; i++) {
    double x_new = x_old + step(gen);
    double p_new = pdf(x_new);
    double p = p_new / p_old;

    if (p < 1) {
      double s = uniform(gen);

      if (s > p) {
        x_new = x_old;
        p_new = p_old;
      }
    }

    x_old = x_new;
    p_old = p_new;

    returner[i] = x_new;
  }

  return returner;
}

double pdf_normal(double x) { return std::exp(-x * x / 2); }

void test_1D_metropolis() {
  using std::cout, std::endl;

  Timer timer;
  auto test_1 = metropolis<double>(pdf_normal, 1000, 0, 0);
  cout << timer.elapsed_ms() << "\n";
  auto test_2 = metropolis(pdf_normal, 1000, 10000, 0);
  cout << timer.elapsed_ms() << "\n";
  auto test_3 = metropolis(pdf_normal, 1000, 0, -4);
  cout << timer.elapsed_ms() << "\n";
  auto test_4 = metropolis(pdf_normal, 1000, 1000, -4);
  cout << timer.elapsed_ms() << endl;
  NumpySaver("build/output/test_normal_distribution.npy")
      << test_1 << test_2 << test_3 << test_4;
}

int main(int argc, char const* argv[]) {
  test_1D_metropolis();
  return 0;
}
