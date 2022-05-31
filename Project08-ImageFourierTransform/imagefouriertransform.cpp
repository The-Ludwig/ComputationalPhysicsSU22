#include <Eigen/Dense>

#include "NumpySaver.hpp"
#include "Timer.hpp"
#include "fft.hpp"
#include "image_fft.hpp"
#include "progressbar.hpp"

using namespace Eigen;
std::tuple<double, double> measure_fft_time(Index input_size) {
  VectorXd vec = VectorXd::Random(input_size);
  VectorXcd dest(input_size);

  return Timer::measure_time(
      [&]() { fft(vec.data(), dest.data(), input_size); });
}

int main(int argc, char const* argv[]) {
#ifndef NDEBUG
  FFT::tests();
  // image_test();
  // compress_test();
#endif

  VectorXd input_sizes_lin = VectorXd::LinSpaced(1030 - 500, 500, 1030);
  VectorXd input_sizes_p2(13);
  // dont judge me
  input_sizes_p2 << 2, 4, 8, 16, 32, 64, 128, 256, 1024, 2048, 4096, 8192,
      16384;

  progressbar bar(input_sizes_lin.size() + input_sizes_p2.size());

  VectorXd mean_lin(input_sizes_lin.size()), std_lin(input_sizes_lin.size());
  for (Index k = 0; k < input_sizes_lin.size(); k++) {
    bar.update();
    auto [mean, std] = measure_fft_time(input_sizes_lin[k]);
    mean_lin[k] = mean;
    std_lin[k] = std;
  }

  VectorXd mean_p2(input_sizes_p2.size()), std_p2(input_sizes_p2.size());
  for (Index k = 0; k < input_sizes_p2.size(); k++) {
    bar.update();
    auto [mean, std] = measure_fft_time(input_sizes_p2[k]);
    mean_p2[k] = mean;
    std_p2[k] = std;
  }

  NumpySaver("build/output/times_lin.npy")
      << input_sizes_lin << mean_lin << std_lin;
  NumpySaver("build/output/times_p2.npy")
      << input_sizes_p2 << mean_p2 << std_p2;
}
