
#include <Eigen/Sparse>
#include <iostream>

#include "NumpySaver.hpp"

using namespace Eigen;

/**
 * @brief Returns the non-zero B-splines at x.
 *
 * @param x The position to evaluate the splines at
 * @param knots Increasing sequence of knots defining the splines. (without
 * ghost points!)
 * @param order The order of the splines this determines the length of the
 * output. indexing convention: lowest order is order=0
 * @return {VectorXd spline_values, idx} the spline values and the index
 */
std::tuple<VectorXd, Index> bsplines(double x, const ArrayXd& knots,
                                     Index order) {
  // TODO: Boundary check on x

  auto n_knots = knots.size();

  auto get_knot = [&knots, &n_knots](int i) {
    if (i < 0)
      return knots[0];
    else if (n_knots <= i)
      return knots[n_knots - 1];
    else
      return knots[i];
  };

  // Find right index
  Index idx = 0;
  while (x > knots(idx + 1)) idx++;

  // check if right ghost points are needed
  if (idx + order > knots.size()) {
  }

  VectorXd iter = VectorXd::Zero(order + 1);
  iter(order) = 1;

  for (Index cur_order = 1; cur_order <= order; cur_order++) {
    for (Index i = idx - cur_order; i < idx; i++) {
      const Index iter_idx = i - idx + order;
      double w1 = 0, w2 = 0;
      // check if both ghost points
      if (i < n_knots && i + cur_order > 0)
        w1 = (x - get_knot(i)) / (get_knot(i + cur_order) - get_knot(i));
      if (i + 1 < n_knots)
        w2 = (get_knot(i + cur_order + 1) - x) /
             (get_knot(i + cur_order + 1) - get_knot(i + 1));

      iter(iter_idx) = w1 * iter(iter_idx) + w2 * iter(iter_idx + 1);
    }
    // this replaces the left 'ghost points'
    double w1 = 0;
    if (idx < knots.size())
      w1 = (x - get_knot(idx)) / (get_knot(idx + cur_order) - get_knot(idx));

    iter(order) = w1 * iter(order);
  }

  return {iter, idx};
}

int main(int argc, char const* argv[]) {
  using std::cout, std::endl;
  ArrayXd knots = ArrayXd::LinSpaced(11, 0, 1);
  cout << "The knots are: " << knots.transpose() << endl;
  double x = .95;
  Index order = 3;
  auto [values, index] = bsplines(x, knots, order);
  cout << "The spline value at " << x << " is: " << values.transpose() << endl;
  cout << "Summing all spline values (should be 1): " << values.sum() << endl;
  return 0;
}
