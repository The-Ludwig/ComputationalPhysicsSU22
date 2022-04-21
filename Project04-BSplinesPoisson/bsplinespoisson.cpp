
#include <fmt/core.h>

#include <Eigen/Sparse>
#include <iostream>
#include <vector>

#include "NumpySaver.hpp"
#include "Timer.hpp"
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

  SparseMatrix<double> matrix(n + 1, n + order - 1);
  matrix.reserve((order + 1) * (n + 1));

  auto [deriv_0, _index] = bsplines::ndx_bsplines(knots[0], knots, order, 2);
  for (int col = 0; col < order - 1; col++)
    matrix.insert(0, col) = deriv_0[col + 1];

  for (int row = 1; row < n; row++) {
    auto [deriv, index] = bsplines::ndx_bsplines(knots[row], knots, order, 2);
    for (auto c : deriv) cout << c << " ";
    cout << "\n";
    for (int col = row - 1; col < row - 1 + order; col++)
      matrix.insert(row, col) = deriv[col - row + 2];
  }

  auto [deriv_last, __index] =
      bsplines::ndx_bsplines(knots[n - 1], knots, order, 1);
  matrix.insert(n, n + order - 2) = deriv_last[order];
  matrix.insert(n, n + order - 3) = deriv_last[order - 1];

  return matrix;
}

int main(int argc, char const* argv[]) {
  cout << "Testing if all spline values at each point sum to 1, no error means "
          "it worked."
       << endl;
  for (double x = 0; x <= 1; x += .001) {
    test_bsplines(x, 3, 11);
  }

  std::vector<double> std_knots(11);
  Map<ArrayXd> knots(std_knots.data(), std_knots.size());
  knots = ArrayXd::LinSpaced(11, 0, 1);

  auto matrix = build_matrix(std_knots);
  cout << matrix.rows() << "x" << matrix.cols() << endl;
  cout << matrix << endl;

  return 0;
}
