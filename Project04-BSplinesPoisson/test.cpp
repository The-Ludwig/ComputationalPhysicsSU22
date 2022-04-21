#include "bspline.hpp"
#include "iostream"

int main() {
  std::vector<double> knots = {0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1};
  // for this we need C++17 structured binding
  int order = 3;
  std::vector<double> xs = {0, .5, 1};

  double x = .5;
  auto [spline_values, index] =
      bspline::ndx_bsplines(x, knots, std::size_t(order), 2);
  for (int i = 0; i <= order; i++)
    std::cout << "Spline second derivative with index "
              << int(index) - int(order - i) << " at " << x << " ="
              << spline_values[i] << std::endl;
  // Expected output:
  // Spline with index 1 (0.5) =0
  // Spline with index 2 (0.5) =0.166667
  // Spline with index 3 (0.5) =0.666667
  // Spline with index 4 (0.5) =0.166667
}
