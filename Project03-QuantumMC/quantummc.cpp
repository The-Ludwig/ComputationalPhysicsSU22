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

double pdf_normal(const Array1d& x) { return std::exp(-x[0] * x[0] / 2); }

void test_1D_metropolis() {
  Array1d argstart;
  argstart << 0;

  std::uniform_real_distribution<> step(-2, 2);
  std::normal_distribution<> step_gauss(0, 1);
  MetropolisAlgorithm<1, std::normal_distribution<> > malg(pdf_normal, argstart,
                                                           step_gauss);

  Timer timer;
  auto test_1 = malg.get_sample(1000);
  cout << timer.elapsed_ms() << "\n";
  auto test_2 = malg.get_sample(1000);
  cout << timer.elapsed_ms() << "\n";
  // auto test_3 = metropolis(pdf_normal, 1000, 0, -4);
  // cout << timer.elapsed_ms() << "\n";
  // auto test_4 = metropolis(pdf_normal, 1000, 1000, -4);
  // cout << timer.elapsed_ms() << endl;
  NumpySaver("build/output/test_normal_distribution.npy") << test_1 << test_2;
}

int main(int argc, char const* argv[]) {
  test_1D_metropolis();
  return 0;
}
