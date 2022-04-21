#pragma once

#include <cassert>
#include <iostream>
#include <tuple>
#include <vector>

namespace bspline {

/**
 * @brief Returns the non-zero B-splines at x. (C++ 17 required, compile with
 * `-std=c++17`)
 *
 * This uses a similar algorithm to DeBoors algorithm
 * (https://en.wikipedia.org/wiki/De_Boor%27s_algorithm) but instead of
 * transforming the recursion and returning the sum of the splines, it returns
 * all the values of the non-zero BSplines at this point.
 *
 *
 * @code
#include <iostream>

#include "bspline.hpp"
...

std::vector<double> knots = {0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1};
// for this we need C++17 structured binding
auto order = 3;
auto x = .5;
auto [spline_values, index] = bspline::bsplines(x, knots, order);
for (int i = 0; i <= order; i++)
  std::cout << "Spline with index " << int(index) - int(order - i) << " ("
            << x << ") =" << spline_values[i] << std::endl;
// Expected output:
// Spline with index 1 (0.5) =0
// Spline with index 2 (0.5) =0.166667
// Spline with index 3 (0.5) =0.666667
// Spline with index 4 (0.5) =0.166667
 * @endcode
 *
 * @tparam T real valued type to work with, e.g. double or float.
 * @param knots Stricly increasing sequence of knots defining the splines.
 * Not including 'ghost points' (starting and trailing repeating points), they
 * are added automaticaly!
 * @param x The position to evaluate the splines at.
 * Must fulfill knots[0] <= x<=knots[knots.size()]
 * @param order The order (often called k) of the splines this determines the
 * length of the output. indexing convention: lowest order is order=0
 * @return {VectorXd spline_values, index} the `order+1` nonzero-spline values
 * and the *last* index of the spline, meaning the splines will have indices in
 * the folowing order: `idx-order, idx-order+1, …, idx-1, idx`.
 * If x is exactly a knot value, `index` is the index to the left of it,
 * except for `x==knots[0]`, then `index==0`.
 */
template <typename T>
std::tuple<std::vector<T>, std::size_t> bsplines(T x,
                                                 const std::vector<T>& knots,
                                                 std::size_t order) {
  int n_knots = knots.size();

  // Check that input fulfills contrains, if in debug mode
#ifndef NDEBUG
  assert((knots[0] <= x) && (x <= knots[n_knots - 1]));
  for (int i = 1; i < n_knots; i++) {
    assert(knots[i - 1] <= knots[i]);
  }
#endif

  // function to get the knot, or automatically return ghost point if our\t of
  // range
  auto get_knot = [&knots, &n_knots](int i) {
    if (i < 0)
      return knots[0];
    else if (n_knots <= i)
      return knots[n_knots - 1];
    else
      return knots[i];
  };

  // Find the knot index where the value is
  int idx = 0;
  while (x > knots[idx + 1]) idx++;

  // initialize the vector with 0, except at the end
  auto iter = std::vector<T>(order + 1, 0);
  iter[order] = 1;

  // fill the vector from idx-cur_order to end (each iteration)
  // the last iteration will contain the non-zero b-spline
  // this algorithm is similar to DeBoor's algorithm
  // but does not calculate the sum, instead, it calculates
  // the single non-zero BSplines
  for (std::size_t cur_order = 1; cur_order <= order; cur_order++) {
    // this replaces the left 'ghost points'
    double w2 = (get_knot(idx + 1) - x) /
                (get_knot(idx + 1) - get_knot(idx - cur_order + 1));
    iter[order - cur_order] = w2 * iter[order - cur_order + 1];

    for (int i = idx - cur_order + 1; i < idx; i++) {
      const int iter_idx = i - idx + order;

      double w1 = (x - get_knot(i)) / (get_knot(i + cur_order) - get_knot(i));
      double w2 = (get_knot(i + cur_order + 1) - x) /
                  (get_knot(i + cur_order + 1) - get_knot(i + 1));

      // iteration
      iter[iter_idx] = w1 * iter[iter_idx] + w2 * iter[iter_idx + 1];
    }

    // this replaces the right 'ghost points'
    double w1 =
        (x - get_knot(idx)) / (get_knot(idx + cur_order) - get_knot(idx));
    // last iteration (iter_idx+1 will always be a ghost point or zero)
    iter[order] = w1 * iter[order];
  }

  // return the bsplines and the index
  return {iter, idx};
}

/**
 * @brief Returns the non-zero second derivatives of B-splines at x. (C++ 17
required, compile with
 * `-std=c++17`)
 *
 * This uses a similar algorithm to DeBoors algorithm
 * (https://en.wikipedia.org/wiki/De_Boor%27s_algorithm) but instead of
 * transforming the recursion and returning the sum of the splines, it returns
 * all the values of the non-zero BSplines at this point.
 *
 *
 * @code
#include <iostream>

#include "bspline.hpp"
...

std::vector<double> knots = {0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1};
// for this we need C++17 structured binding
auto order = 3;
auto x = .5;
auto [spline_values, index] = bspline::bsplines(x, knots, order);
for (int i = 0; i <= order; i++)
  std::cout << "Spline with index " << int(index) - int(order - i) << " ("
            << x << ") =" << spline_values[i] << std::endl;
// Expected output:
// Spline with index 1 (0.5) =0
// Spline with index 2 (0.5) =0.166667
// Spline with index 3 (0.5) =0.666667
// Spline with index 4 (0.5) =0.166667
 * @endcode
 *
 * @tparam T real valued type to work with, e.g. double or float.
 * @param knots Stricly increasing sequence of knots defining the splines.
 * Not including 'ghost points' (starting and trailing repeating points), they
 * are added automaticaly!
 * @param x The position to evaluate the splines at.
 * Must fulfill knots[0] <= x<=knots[knots.size()]
 * @param order The order (often called k) of the splines this determines the
 * length of the output. indexing convention: lowest order is order=0
 * @return {VectorXd spline_values, index} the `order+1` nonzero-spline values
 * and the *last* index of the spline, meaning the splines will have indices in
 * the folowing order: `idx-order, idx-order+1, …, idx-1, idx`.
 * If x is exactly a knot value, `index` is the index to the left of it,
 * except for `x==knots[0]`, then `index==0`.
 */
template <typename T>
std::tuple<std::vector<T>, std::size_t> dxbsplines(T x,
                                                   const std::vector<T>& knots,
                                                   std::size_t order) {
  int n_knots = knots.size();

  // Check that input fulfills contrains, if in debug mode
#ifndef NDEBUG
  assert((knots[0] <= x) && (x <= knots[n_knots - 1]));
  for (int i = 1; i < n_knots; i++) {
    assert(knots[i - 1] <= knots[i]);
  }
#endif

  // function to get the knot, or automatically return ghost point if our\t of
  // range
  auto get_knot = [&knots, &n_knots](int i) {
    if (i < 0)
      return knots[0];
    else if (n_knots <= i)
      return knots[n_knots - 1];
    else
      return knots[i];
  };

  // Find the knot index where the value is
  int idx = 0;
  while (x > knots[idx + 1] && idx < n_knots - 1) idx++;

  // initialize the vector with 0, except at the end
  auto iter = std::vector<T>(order + 1, 0);
  iter[order] = 1;

  // fill the vector from idx-cur_order to end (each iteration)
  // the last iteration will contain the non-zero b-spline
  // this algorithm is similar to DeBoor's algorithm
  // but does not calculate the sum, instead, it calculates
  // the single non-zero BSplines
  for (std::size_t cur_order = 1; cur_order <= order; cur_order++) {
    for (int i = idx - int(cur_order); i < idx; i++) {
      const int iter_idx = i - idx + int(order);
      double w1 = 0, w2 = 0;

      // check if both ghost points, to not devide by zero!
      if (i < n_knots && i + cur_order > 0)
        w1 = (x - get_knot(i)) / (get_knot(i + cur_order) - get_knot(i));
      if (i + 1 < n_knots)
        w2 = (get_knot(i + cur_order + 1) - x) /
             (get_knot(i + cur_order + 1) - get_knot(i + 1));

      // iteration
      iter[iter_idx] = w1 * iter[iter_idx] + w2 * iter[iter_idx + 1];
    }

    // this replaces the left 'ghost points'
    double w1 =
        (x - get_knot(idx)) / (get_knot(idx + cur_order) - get_knot(idx));
    // last iteration (iter_idx+1 will always be a ghost point or zero)
    iter[order] = w1 * iter[order];
  }

  // return the bsplines and the index
  return {iter, idx};
}

}  // namespace bspline