#include <yaml-cpp/yaml.h>

#include <iostream>

#include "NumpySaver.hpp"

using std::size_t, std::vector;

template <class T>
void print_vec(std::vector<T>& vec, const char* name) {
  std::cout << name << " = ";
  for (auto& v : vec) std::cout << v << ", ";
  std::cout << std::endl;
}

void solve(vector<double>& hs, double incoming_radiation, double alpha_vis,
           double alpha_ir, std::size_t iters = 1000) {
  using namespace std;
  size_t N = hs.size();

  double kelvin_0 = 273.15;

  double rho_b = 1.224;        // kg/m3
  double t_b = 20 + kelvin_0;  // K
  double h_b = hs[0];
  double g0 = 9.81;
  double M = 0.02896;
  double r = 8.3145;

  double sigma_sb = 5.67e-8;  // W/K^4/m^2

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

  // initialize the iteration vectors
  vector<double> t_in_ir(N, 0);
  vector<double> t_out(N, 0);
  vector<double> e(N, 0);

  t_in_ir[N - 1] = 0;
  t_out[N - 1] = 0;
  e[N - 1] = 0;
  t_out[0] = 0;

  // iterate
  for (size_t k = 0; k < iters; k++) {
    for (size_t i = N - 2; i > 0; i--) {
      double h = hs[i + 1] - hs[i];
      t_in_ir[i] =
          (t_in_ir[i + 1] + e[i + 1] / 2.) * exp(-sigma_ir * rho[i] * h);
      t_out[i] = (t_out[i - 1] + e[i - 1] / 2.) * exp(-sigma_ir * rho[i] * h);
      e[i] = (.5 * (e[i - 1] + e[i + 1]) + t_out[i - 1] + t_in_ir[i + 1]) *
                 (1. - exp(-sigma_ir * rho[i] * h)) +
             t_in_v[i + 1] * (1. - exp(-sigma_vis * rho[i] * h));
    }
    double h = hs[1] - hs[0];
    t_in_ir[0] = (t_in_ir[1] + e[1] / 2.) * exp(-sigma_ir * rho[0] * h);
    e[0] = (.5 * e[1] + t_in_ir[1]) * (1. - exp(-sigma_ir * rho[0] * h)) +
           t_in_v[1] * (1. - exp(-sigma_vis * rho[0] * h));

    cout << "T = " << sqrt(sqrt(e[0] / sigma_sb)) - kelvin_0 << " C" << endl;
  }
}

int main() {
  using namespace std;

  size_t N = 10;
  vector<double> hs(N);
  double h0 = 0;         // m
  double hmax = 100000;  // m
  for (size_t i = 0; i < hs.size(); i++)
    hs[i] = h0 + (hmax - h0) / (double(N) - 1.) * double(i);

  // print_vec(hs, "h = ");

  solve(hs, 344. * (1. - .3), 1e-4, 1e-1);
  return 0;
}
