#include <fmt/core.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <tuple>
#include <vector>

#include "NumpySaver.hpp"

using std::size_t, std::vector;

template <class T>
void print_vec(std::vector<T>& vec, const char* name) {
  std::cout << name << " = ";
  for (auto& v : vec) std::cout << v << ", ";
  std::cout << std::endl;
}

double dif(std::vector<double>& v1, std::vector<double>& v2) {
  using namespace std;
  assert(v1.size() == v2.size());
  double sum = 0;
  for (size_t i = 0; i < v1.size(); i++) {
    double t = v1[i] - v2[i];
    sum += t * t;
  }
  return sqrt(sum);
}

std::tuple<double, std::size_t> solve(vector<double>& hs,
                                      double incoming_radiation,
                                      double alpha_vis, double alpha_ir,
                                      std::optional<std::string> basename = {},
                                      double tol = 1e-3) {
  using namespace std;
  size_t N = hs.size();

  constexpr double kelvin_0 = 273.15;
  constexpr double rho_b = 1.224;        // kg/m3
  constexpr double t_b = 20 + kelvin_0;  // K
  double h_b = hs[0];
  constexpr double g0 = 9.81;
  constexpr double M = 0.02896;
  constexpr double r = 8.3145;

  constexpr double sigma_sb = 5.67e-8;  // W/K^4/m^2

  double sigma_vis = alpha_vis / rho_b;
  double sigma_ir = alpha_ir / rho_b;

  // calculate the rho
  vector<double> rho(N);
  for (size_t i = 0; i < N; i++)
    rho[i] = rho_b * exp(-g0 * M * (hs[i] - h_b) / (r * t_b));

  // print_vec(rho, "ρ");

  // calculate the visible incoming energy
  vector<double> t_in_v(N);
  t_in_v[N - 1] = incoming_radiation;
  for (int i = N - 2; i >= 0; i--)
    t_in_v[i] = t_in_v[i + 1] * exp(-(hs[i + 1] - hs[i]) * sigma_vis);

  // print_vec(t_in_v, "T_in_v");

  // initialize the iteration vectors
  vector<double> t_in_ir(N, 0);
  vector<double> t_out(N, 0);
  vector<double> e(N, 0);

  t_in_ir[N - 1] = 0;
  t_out[N - 1] = 0;
  e[N - 1] = 0;
  t_out[0] = 0;
  t_in_ir[0] = 0;

  vector<double> old_e(e.size(), numeric_limits<double>::max());

  size_t iters = 0;
  // iterate
  while (dif(old_e, e) > tol) {
    iters++;
    old_e = e;

    // down
    for (size_t i = N - 2; 1 < i; i--) {
      double h = hs[i + 1] - hs[i];
      t_in_ir[i] =
          (t_in_ir[i + 1] + e[i + 1] / 2.) * exp(-sigma_ir * rho[i] * h);
      t_out[i] = (t_out[i - 1] + e[i - 1] / 2.) * exp(-sigma_ir * rho[i] * h);
      e[i] = (.5 * (e[i - 1] + e[i + 1]) + t_out[i - 1] + t_in_ir[i + 1]) *
                 (1. - exp(-sigma_ir * rho[i] * h)) +
             t_in_v[i + 1] * (1. - exp(-sigma_vis * rho[i] * h));
    }
    double h = hs[2] - hs[1];
    t_out[1] = (t_out[0] + e[0]) * exp(-sigma_ir * rho[1] * h);
    t_in_ir[1] = (t_in_ir[2] + e[2] / 2.) * exp(-sigma_ir * rho[1] * h);

    e[1] = (e[0] + .5 * e[2] + t_out[0] + t_in_ir[2]) *
               (1. - exp(-sigma_ir * rho[1] * h)) +
           t_in_v[2] * (1. - exp(-sigma_vis * rho[1] * h));

    e[0] = (.5 * e[1] + t_in_ir[1]) + t_in_v[1];
  }

  vector<double> ts;
  ts.reserve(e.size());
  for (auto& e_now : e) ts.push_back(sqrt(sqrt(e_now / sigma_sb)) - kelvin_0);

  if (basename.has_value())
    NumpySaver(fmt::format("build/output/{}.npy", basename.value()))
        << hs << rho << e << t_in_v << ts;

  double T = sqrt(sqrt(e[0] / sigma_sb)) - kelvin_0;
  cout << "T (" << iters << " iterations) = " << T << " °C" << endl;

  cout << "Outgoing radiation = " << e[N - 2] / 2. + t_out[N - 2] << endl;

  return {T, iters};
}

std::tuple<double, std::size_t> solve_lin_ground(
    std::size_t N, double incoming_radiation, double alpha_vis, double alpha_ir,
    std::optional<std::string> basename = {}, double tol = 1e-3) {
  vector<double> hs(N);
  double h0 = 0;        // m
  double hmax = 80000;  // m
  hs[0] = h0;
  hs[1] = h0;  // "ground"
  for (size_t i = 1; i < hs.size() - 1; i++)
    hs[i + 1] = h0 + (hmax - h0) / (double(N) - 2.) * double(i);

  return solve(hs, incoming_radiation, alpha_vis, alpha_ir, basename, tol);
}

std::tuple<double, std::size_t> solve_lin(
    std::size_t N, double incoming_radiation, double alpha_vis, double alpha_ir,
    std::optional<std::string> basename = {}, double tol = 1e-3) {
  vector<double> hs(N);
  double h0 = 0;        // m
  double hmax = 80000;  // m
  hs[0] = h0;
  for (size_t i = 1; i < hs.size(); i++)
    hs[i] = h0 + (hmax - h0) / (double(N) - 1.) * double(i);

  return solve(hs, incoming_radiation, alpha_vis, alpha_ir, basename, tol);
}

int main() {
  using namespace std;

  string name("test");
  solve_lin_ground(1000, 344. * (1. - .3), 5e-5, 1e-2, name);

  {
    vector<double> ns = {5, 10, 100, 1000, 10000, 100000};
    vector<double> iters, temps;
    iters.reserve(ns.size());
    temps.reserve(ns.size());
    for (auto& n : ns) {
      string name = fmt::format("N_{}", n);
      auto [temp, iter] =
          solve_lin_ground(n, 344. * (1. - .3), 5e-5, 1e-2, name);
      temps.push_back(temp);
      iters.push_back(iter);
    }
    NumpySaver("build/output/ns.npy") << ns << iters << temps;
  }

  {
    vector<double> sigma_viss = {1e-5, 2e-5, 3e-5, 4e-5, 5e-5,
                                 6e-5, 7e-5, 8e-5, 9e-5, 1e-4};
    vector<double> iters, temps;
    iters.reserve(sigma_viss.size());
    temps.reserve(sigma_viss.size());
    for (auto& sigma_vis : sigma_viss) {
      string name = fmt::format("sigma_{}", sigma_vis);
      auto [temp, iter] =
          solve_lin_ground(1000, 344. * (1. - .3), sigma_vis, 1e-2, name);
      temps.push_back(temp);
      iters.push_back(iter);
    }
    NumpySaver("build/output/sigma_viss.npy") << sigma_viss << iters << temps;
  }

  {
    vector<double> sigma_irs = {1e-3, 3e-3, 5e-3, 1e-2, 5e-2, 10e-2};
    vector<double> iters, temps;
    iters.reserve(sigma_irs.size());
    temps.reserve(sigma_irs.size());
    for (auto& sigma_ir : sigma_irs) {
      string name = fmt::format("sigma_ir_{}", sigma_ir);
      auto [temp, iter] =
          solve_lin_ground(1000, 344. * (1. - .3), 5e-5, sigma_ir, name);
      temps.push_back(temp);
      iters.push_back(iter);
    }
    NumpySaver("build/output/sigma_irs.npy") << sigma_irs << iters << temps;
  }
}
