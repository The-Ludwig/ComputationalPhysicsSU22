
#include <Eigen/Sparse>
#include <iostream>

#include "NumpySaver.hpp"

using namespace Eigen;

/**
 * @brief Returns the non-zero B-splines at x.
 *
 * @param x The position to evaluate the splines at
 * @param knots The knots defining the splines. (without ghost points!)
 * @param order The order of the splines this determines the length of the
 * output. indexing convention: lowest order is order=0
 * @return SparseVector<double> with `order+1` nonzero elements
 */
SparseVector<double> bsplines(double x, const VectorXd& knots, Index order) {
  // TODO: Boundary check on x

  // Find right index
  Index idx = 0;
  while (x > knots(idx + 1)) idx++;

  // check if right ghost points are needed
  if (idx + order > knots.size()) {
  }

  VectorXd iter = VectorXd::Zero(order + 1);
  iter(order) = 1;

  for (Index cur_order = 1; cur_order <= order; cur_order++) {
    for (Index i = order - cur_order; i < order; i++) {
      const Index knots_idx = idx + i - order;

      // check if ghost point
      if (knots_idx >= 0 && knots_idx + cur_order < knots.size()) {
        double w1 = (x - knots[knots_idx]) /
                    (knots[knots_idx + cur_order] - knots[knots_idx]);
        iter(i) = w1 * iter(i);

        if (knots_idx + cur_order + 1 < knots.size()) {
          double w2 = (knots[knots_idx + cur_order + 1] - x) /
                      (knots[knots_idx + cur_order + 1] - knots[knots_idx + 1]);
          iter(i) += w2 * iter(i + 1);
        }
      }
    }
    // this replaces the left 'ghost points'
    double w1 = (x - knots[idx]) / (knots[idx + cur_order] - knots[idx]);
    iter(order) = w1 * iter(order);
  }

  SparseVector<double> return_vector(knots.size());
  return_vector.reserve(order + 1);
  for (Index i = idx - order; i < idx; i++)
    return_vector.insert(i) = iter(i - idx + order);

  return return_vector;
}

int main(int argc, char const* argv[]) {
  Eigen::ArrayXd array;
  std::string name;

  return 0;
}
