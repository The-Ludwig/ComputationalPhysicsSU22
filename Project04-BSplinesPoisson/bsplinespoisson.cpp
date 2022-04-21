
#include <fmt/core.h>

#include <Eigen/Sparse>
#include <iostream>

#include "NumpySaver.hpp"
#include "Timer.hpp"

using namespace Eigen;

// void test_bsplines(double x = .1, Index order = 3, Index nknots = 11) {
//   using std::cout, std::endl, std::cerr;

//   ArrayXd knots = ArrayXd::LinSpaced(nknots, 0, 1);
//   Timer t;
//   auto [values, index] = bsplines(x, knots, order);
//   cout << fmt::format(
//       "Testing BSpline algorithm with order={}, {} knots took {:.3f}ms\n",
//       (int)order, (int)nknots, t.elapsed_ms());

//   auto sum = values.sum();
//   if (std::abs(sum - 1) > 1e-10) {
//     cerr << "The knots are: " << knots.transpose() << endl;
//     cerr << "The spline value at " << x << " is: " << values.transpose()
//          << endl;
//     cerr << "Summing all spline values (should be 1): " << sum << endl;
//   }
// }

int main(int argc, char const* argv[]) {
  // for (int i = 11; i < 5000; i += 100) {
  //   test_bsplines(.5, 3, i);
  // }

  // ArrayXd knots = ArrayXd::LinSpaced(11, 0, 1);
  // auto [values, index] = bsplines(.1, knots, 3);
  // std::cout << values;

#ifdef NDEBUG
  std::cout << "This is not a debug session" << std::endl;
#else
  std::cout << "This is a debug session" << std::endl;
#endif

  return 0;
}
