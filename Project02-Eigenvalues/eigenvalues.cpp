#include <yaml-cpp/yaml.h>

#include <Eigen/Dense>
#include <Eigen/OrderingMethods>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <iostream>
#include <vector>

#include "NumpySaver.hpp"

using namespace Eigen;
using namespace std;

typedef double type;

const char* divider = "###################################\n";

// copied from https://gist.github.com/lorenzoriano/5414671
template <typename T>
std::vector<T> linspace(T start, T end, size_t N) {
  T h = (end - start) / static_cast<T>(N - 1);
  std::vector<T> xs(N);
  typename std::vector<T>::iterator x;
  T val;
  for (x = xs.begin(), val = start; x != xs.end(); ++x, val += h) *x = val;
  return xs;
}

template <typename T>
vector<T> get_values(YAML::Node& node) {
  auto type = node.Type();
  // single value
  if (type == YAML::NodeType::Scalar) {
    vector<T> ret = {node.as<T>()};
    return ret;
  } else if (type == YAML::NodeType::Sequence) {
    return node.as<vector<T>>();
  } else if (type == YAML::NodeType::Map) {
    auto start = node["start"].as<T>();
    auto end = node["end"].as<T>();
    auto N = node["N"].as<size_t>();
    return linspace<T>(start, end, N);
  } else {
    throw std::runtime_error(
        "Range node must either be list, scalar or a map defining "
        "'start','end' and 'N'.");
  }
}

tuple<VectorXd, type, unsigned int> eigenvector_inverse_power_iteration(
    SparseMatrix<type>& matrix, VectorXd& start_vector, double shift = 0,
    double tol = 1e-4, unsigned int max_iterations = 10000) {
  // apply shift
  // cout << "Matrix=\n" << matrix;
  SparseMatrix<type> mat(matrix.rows(), matrix.cols());
  mat.setIdentity();
  mat = matrix - shift * mat;
  // cout << "Shifted Matrix\n" << mat;

  // first do the (sparse) LU decomposition
  Eigen::SparseLU<Eigen::SparseMatrix<type>, Eigen::COLAMDOrdering<int>> solver;
  solver.analyzePattern(mat);
  solver.factorize(mat);

  // now the iteration
  VectorXd y_now(start_vector);
  VectorXd y_last;

  // not the best init value, but should work ok
  type lambda_now = numeric_limits<type>::max();
  type lambda_last;

  // keep track of iterations
  unsigned int i = 0;

  do {
    // normalize vector
    y_now.normalize();
    y_last = y_now;

    // iteration
    y_now = solver.solve(y_last);

    // calculate eigenvector
    lambda_last = lambda_now;
    lambda_now = 1 / (y_last.dot(y_now)) + shift;

    if (i++ >= max_iterations) break;
  } while (abs(lambda_last - lambda_now) > tol);

  return {y_now, lambda_now, i};
}

void test_inverse_power_iteration(YAML::Node& test_configs) {
  auto size = test_configs.size();
  for (std::size_t i = 0; i < size; i++) {
    auto conf = test_configs[i];

    auto N = conf["dim"].as<int>();
    auto max_iteartions = conf["max_iter"].as<unsigned int>();
    auto tol = conf["tol"].as<double>();

    MatrixXd matrix = MatrixXd::Random(N, N);
    matrix = matrix.selfadjointView<Eigen::Lower>();
    VectorXd start_vector = VectorXd::Random(N);
    SparseMatrix<type> sparse_view = matrix.sparseView();

    cout << divider << "# Test Power Iteration " << i + 1 << "/" << size << endl
         << divider;
    cout << "Matrix=\n" << matrix;
    cout << "\nStart vector=\n" << start_vector;
    cout << "\nTolerance=" << tol << endl << endl;

    YAML::Node shift_node = conf["shift"];
    auto shifts = get_values<double>(shift_node);

    // Eigen solution for comparison:
    EigenSolver<MatrixXd> es(matrix);
    cout << "Eigen-Eigenvalues:" << endl << es.eigenvalues() << endl;
    cout << "Eigen-Eigenvectors:" << endl << es.eigenvectors() << endl << endl;

    for (auto& shift : shifts) {
      cout << divider << "Shift=" << shift << endl;

      auto [eigenvector, eigenvalue, iterations] =
          eigenvector_inverse_power_iteration(sparse_view, start_vector, shift,
                                              tol, max_iteartions);
      cout << "Iterations needed to reach tolerance: " << iterations << endl;
      cout << "Inverse power Iteration eigenvector is:\n" << eigenvector;
      cout << "\nInverse power Iteration eigenvalue is:\n" << eigenvalue;
      if (iterations >= max_iteartions)
        cerr
            << "\nWarning! Exceeded maximum number of iterations! Solution did "
               "not converge!"
            << endl;
      cout << endl;
    }
    cout << endl;
  }
}

void run_tests() {
  cout << divider << divider << "## RUNNING TESTS\n"
       << divider << divider << endl
       << endl;

  YAML::Node config = YAML::LoadFile("test.yaml");

  auto config_power_iter = config["inverse_iteration"];
  test_inverse_power_iteration(config_power_iter);
}

int main(int argc, char const* argv[]) {
  run_tests();
  return 0;
}
