#include <yaml-cpp/yaml.h>

#include <Eigen/Dense>
#include <Eigen/OrderingMethods>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <chrono>
#include <iostream>
#include <vector>

#include "NumpySaver.hpp"

using namespace Eigen;
using namespace std;

typedef double type;

const char* divider = "#####################################################";

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

SparseMatrix<type> get_5_point_schroedinger_matrix(
    double xmin, double xmax, unsigned int N,
    function<double(double)> potential) {
  // cout << "Hello from the Beginning of building the matrix" << endl;
  SparseMatrix<type> mat(N, N);
  mat.reserve(5 * N - 6);

  double h = (xmax - xmin) / N;
  double prefac = 1. / 12. / h / h;

  double x = xmin;
  for (unsigned int i = 0; i < N; i++, x += h) {
    if (i >= 2) mat.insert(i, i - 2) = prefac;
    if (i >= 1) mat.insert(i, i - 1) = -prefac * 16;
    mat.insert(i, i) = prefac * 30 + potential(x);
    if (i < N - 1) mat.insert(i, i + 1) = -prefac * 16;
    if (i < N - 2) mat.insert(i, i + 2) = prefac;
  }
  // cout << "Hello from the End of building the matrix" << endl;

  return mat;
}

double harmonic_oszi_potential(double x) { return x * x; }
double disturbed_harmonic_oszi_potential(double x, double c1, double c2) {
  return x * x + c1 * exp(-x * x * c2);
}

function<double(double)> get_potential(YAML::Node potential_conf) {
  auto pot_name = potential_conf["name"].as<string>();

  if (pot_name == "harmonic") {
    return harmonic_oszi_potential;
  } else if (pot_name == "disturbed") {
    double c1 = potential_conf["c1"].as<double>();
    double c2 = potential_conf["c2"].as<double>();
    return bind(disturbed_harmonic_oszi_potential, std::placeholders::_1, c1,
                c2);
  } else {
    throw runtime_error(
        "Potential needs to be either 'harmonic' or 'disturbed'!");
  }
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

  } while ((abs(lambda_last - lambda_now) > tol) && i++ < max_iterations);

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

void test_harmonic(YAML::Node& test_configs) {
  auto size = test_configs.size();
  for (std::size_t i = 0; i < size; i++) {
    cout << divider << "\n# Test Harmonic Oscillator " << i + 1 << "/" << size
         << endl
         << divider << endl
         << endl;
    auto conf = test_configs[i];

    auto N = conf["N"].as<size_t>();
    auto max_iteartions = conf["max_iter"].as<unsigned int>();
    auto tol = conf["tol"].as<double>();
    auto xmin = conf["xmin"].as<double>();
    auto xmax = conf["xmax"].as<double>();

    function<double(double)> pot = get_potential(conf["potential"]);
    auto matrix = get_5_point_schroedinger_matrix(xmin, xmax, N, pot);

    bool has_name = false;
    string file_name;
    if (conf["name"]) {
      has_name = true;
      file_name = "build/output/" + conf["name"].as<string>() + ".npy";
      cout << "Name to save to: " << file_name << endl;
    }

    NumpySaver saver(file_name);

    if (has_name) {
      VectorXd x_space(N);
      VectorXd pot_space(N);
      double h = (xmax - xmin) / N;

      double x = xmin;
      for (size_t i = 0; i < N; i++, x += h) {
        x_space[i] = x;
        pot_space[i] = pot(x);
      }

      saver << "x" << x_space << "pot" << pot_space;
    }

    VectorXd start_vector = VectorXd::Ones(N);

    YAML::Node shift_node = conf["shift"];
    auto shifts = get_values<double>(shift_node);

    cout << "Tolerance=" << tol << endl;

    // Eigen solution for comparison:
    SelfAdjointEigenSolver<SparseMatrix<type>> es(matrix);
    if (N < 10) {
      cout << "Eigen-Eigenvalues:" << endl << es.eigenvalues() << endl;
      cout << "Eigen-Eigenvectors:" << endl
           << es.eigenvectors() << endl
           << endl;
      cout << "Matrix=\n" << matrix << endl;
      cout << "Start vector=\n" << start_vector << endl;
    }

    for (auto& shift : shifts) {
      cout << divider << "Shift=" << shift << endl;

      auto [eigenvector, eigenvalue, iterations] =
          eigenvector_inverse_power_iteration(matrix, start_vector, shift, tol,
                                              max_iteartions);

      if (has_name) saver << "ev" << eigenvector;

      cout << "Iterations=" << iterations << endl;

      if (N < 10)
        cout << "Inverse power Iteration eigenvector = \n"
             << eigenvector << endl;

      cout << "Inverse power Iteration eigenvalue = " << eigenvalue << endl;

      if (iterations >= max_iteartions)
        cerr << "Warning! Exceeded maximum number of iterations! Solution did "
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

  auto config_harmonic = config["harmonic"];
  test_harmonic(config_harmonic);
}

void run_program() {
  YAML::Node configs = YAML::LoadFile("config.yaml");
  auto size = configs.size();

  for (std::size_t i = 0; i < size; i++) {
    cout << divider << divider << "\n# Running Simulation  " << i + 1 << "/"
         << size << endl
         << divider << divider << endl;
    auto conf = configs[i];

    auto N = conf["N"].as<size_t>();
    auto max_iteartions = conf["max_iter"].as<unsigned int>();
    auto tol = conf["tol"].as<double>();
    auto xmin = conf["xmin"].as<double>();
    auto xmax = conf["xmax"].as<double>();

    function<double(double)> pot = get_potential(conf["potential"]);
    auto matrix = get_5_point_schroedinger_matrix(xmin, xmax, N, pot);
    VectorXd start_vector = VectorXd::Ones(N);

    bool has_name = false;
    string file_name;
    if (conf["name"]) {
      has_name = true;
      file_name = "build/output/" + conf["name"].as<string>() + ".npy";
      cout << "Save to:\t" << file_name << endl;
    }

    NumpySaver saver(file_name);

    if (has_name) {
      VectorXd x_space(N);
      VectorXd pot_space(N);
      double h = (xmax - xmin) / N;

      double x = xmin;
      for (size_t i = 0; i < N; i++, x += h) {
        x_space[i] = x;
        pot_space[i] = pot(x);
      }

      saver << "x" << x_space << "pot" << pot_space;
    }

    YAML::Node shift_node = conf["shift"];
    auto shifts = get_values<double>(shift_node);

    cout << "N=\t\t" << N << endl;
    cout << "Tolerance=\t" << tol << endl;
    cout << "xmin=\t\t" << xmin << endl;
    cout << "xmax=\t\t" << xmax << endl;
    cout << "function=\t" << conf["potential"]["name"] << endl << endl;

    // to time the methods
    chrono::steady_clock::time_point start;
    chrono::duration<double> elapsed_seconds;

    // Eigen solution for comparison:
    if (N <= 1000) {
      start = chrono::steady_clock::now();

      SelfAdjointEigenSolver<SparseMatrix<type>> es(matrix);
      auto ev_ref = es.eigenvalues();

      elapsed_seconds = chrono::steady_clock::now() - start;
      cout << "Eigen::SelfAdjointEigenSolver<SparseMatrix> took "
           << elapsed_seconds.count() * 1e3 << "ms" << endl;

      cout << "eigenvalues=" << ev_ref(seq(0, shifts.size())).transpose()
           << endl
           << endl;

      auto evec_ref = es.eigenvectors();

      if (has_name)
        NumpySaver("build/output/" + conf["name"].as<string>() +
                   "_reference.npy")
            << evec_ref;

      if (N < 10) {
        cout << "Eigen-Eigenvalues:" << endl << ev_ref << endl;
        cout << "Eigen-Eigenvectors:" << endl << evec_ref << endl << endl;
        cout << "Matrix=\n" << matrix << endl;
        cout << "Start vector=\n" << start_vector << endl;
      }
    } else {
      cout << "Not going to test SelfAdjointEigenSolver, since it would take "
              "too long."
           << endl;
    }

    start = std::chrono::steady_clock::now();
    for (auto& shift : shifts) {
      cout << divider << "\nShift=\t\t" << shift << endl;

      auto [eigenvector, eigenvalue, iterations] =
          eigenvector_inverse_power_iteration(matrix, start_vector, shift, tol,
                                              max_iteartions);

      elapsed_seconds = chrono::steady_clock::now() - start;

      if (has_name) saver << "ev" << eigenvector;

      cout << "elapsed time=\t" << elapsed_seconds.count() * 1e3 << "ms"
           << endl;

      cout << "Iterations=\t" << iterations << endl;
      if (N < 10)
        cout << "Inverse power Iteration eigenvector = \n"
             << eigenvector << endl;

      cout << "eigenvalue=\t" << eigenvalue << endl;

      if (iterations >= max_iteartions)
        cerr << "Warning! Exceeded maximum number of iterations! Solution did "
                "not converge!"
             << endl;
    }
    cout << endl;
  }
}

int main(int argc, char const* argv[]) {
  // uncomment this to run the test
  // run_tests();
  run_program();
  return 0;
}
