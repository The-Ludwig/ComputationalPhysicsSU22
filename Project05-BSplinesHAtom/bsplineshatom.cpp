#include <yaml-cpp/yaml.h>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <boost/math/quadrature/gauss.hpp>
#include <iostream>
#include <vector>

#include "NumpySaver.hpp"
#include "bsplines.hpp"

using namespace Eigen;
using std::cout, std::endl;

template <int num_points>
double matrix_B_element(int i, int j, std::vector<double>& knots, Index order) {
  using boost::math::quadrature::gauss;
  int imax = std::max(i, j);
  int imin = std::min(i, j);

  int dif = imax - imin;

  if (dif >= order) return 0;

  Index n = knots.size();

  // alternatively do it between each knot
  // double sum = 0;
  // for (int i = std::max(0, imax); i < std::min(imin + order, n - 1); i++) {
  //   cout << i << "-" << i + 1 << "=" << knots[i + 1] - knots[i] << endl;
  // auto f = [&](const double& t) {
  //     auto [values, idx] = bsplines::bsplines(t, knots, order);
  //     return values[imax - ((int)idx - (int)order)] *
  //            values[imin - ((int)idx - (int)order)];
  //   };
  //   sum += gauss<double, num_points>::integrate(
  //       f, knots[i], knots[i+1]);
  // }

  auto f = [&](const double& t) {
    auto [values, idx] = bsplines::bsplines(t, knots, order);
    return values[imax - ((int)idx - (int)order)] *
           values[imin - ((int)idx - (int)order)];
  };
  double sum = gauss<double, num_points>::integrate(
      f, knots[std::max(0, imax)], knots[std::min(imin + order, n - 1)]);

  return sum;
}

template <int num_points>
double matrix_H_element(int i, int j, std::vector<double>& knots, Index order,
                        int l, int z) {
  using boost::math::quadrature::gauss;
  int imax = std::max(i, j);
  int imin = std::min(i, j);

  int dif = imax - imin;

  if (dif > order) return 0;

  Index n = knots.size();

  auto f = [&](const double& t) {
    auto [values, idx] = bsplines::bsplines(t, knots, order);
    auto [values_dd, idx_dd] = bsplines::ndx_bsplines(t, knots, order, 2);
    double h_val =
        -values_dd[j - (idx - order)] / 2. +
        (l * (l + 1.) / (t * t * 2.) - z / t) * values[j - (idx - order)];
    return values[i - (idx - order)] * h_val;
    ;
  };

  double sum = gauss<double, num_points>::integrate(
      f, knots[std::max(0, imax)], knots[std::min(imin + order, n - 1)]);

  return sum;
}

template <int num_points>
MatrixXd matrix_B(std::vector<double>& knots, Index order) {
  Index n = knots.size();
  MatrixXd mat(n - 3 + order, n - 3 + order);

  for (int row = 0; row < mat.rows(); row++) {
    for (int col = 0; col < mat.cols(); col++) {
      mat(row, col) = matrix_B_element<num_points>(
          row - order + 1, col - order + 1, knots, order);
    }
  }

  return mat;
}

template <int num_points>
MatrixXd matrix_H(std::vector<double>& knots, Index order, int l, int z) {
  Index n = knots.size();
  MatrixXd mat(n - 3 + order, n - 3 + order);

  for (int row = 0; row < mat.rows(); row++) {
    for (int col = 0; col < mat.cols(); col++) {
      mat(row, col) = matrix_H_element<num_points>(
          row - order + 1, col - order + 1, knots, order, l, z);
    }
  }

  return mat;
}

int main() {
  int n_knots = 15;
  std::vector<double> knots(n_knots);
  Map<ArrayXd> eigen_knots(knots.data(), knots.size());
  eigen_knots = ArrayXd::LinSpaced(n_knots, 0, 10);

  cout << "knots = " << eigen_knots.transpose() << endl;

  auto B = matrix_B<3000>(knots, 3);
  auto H = matrix_H<3000>(knots, 3, 0, 1);

  cout << "B=\n" << B << "\n\n";
  cout << "H=\n" << H << endl;

  GeneralizedEigenSolver<MatrixXd> ges;
  ges.compute(H, B);
  cout << "Eigenvalues: " << ges.eigenvalues().transpose() << endl;

  return 0;
}
