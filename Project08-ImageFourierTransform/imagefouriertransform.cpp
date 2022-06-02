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
      [&]() { fft(vec.data(), dest.data(), input_size); }, .01);
}

int main(int argc, char const* argv[]) {
#ifndef NDEBUG
  FFT::tests();
  // image_test();
  compress_test();
#endif

  VectorXd input_sizes_lin = VectorXd::LinSpaced(1030 - 500 + 1, 500, 1030);
  VectorXd input_sizes_p2(14);
  // dont judge me
  input_sizes_p2 << 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192,
      16384;

  progressbar bar(input_sizes_lin.size() + input_sizes_p2.size());

  VectorXd mean_lin(input_sizes_lin.size()), std_lin(input_sizes_lin.size());
  for (Index k = 0; k < input_sizes_lin.size(); k++) {
    bar.update();
    auto [mean, std] = measure_fft_time(input_sizes_lin[k] + 1e-5);
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

  VectorXd lost(6);
  VectorXcd lost_fft(lost.size());
  VectorXd lost_dct(lost.size());
  lost << 4, 8, 15, 16, 23, 42;
  fft(lost.data(), lost_fft.data(), lost.size());
  dct(lost.data(), lost_dct.data(), lost.size());
  NumpySaver("build/output/lost.npy")
      << lost << lost_fft.real() << lost_fft.imag() << lost_dct;

  VectorXd plot(1000);
  VectorXcd plot_fft(plot.size());
  for (Index k = 0; k < plot.size() / 3; k++) {
    plot[k] = 0;
  }
  for (Index k = plot.size() / 3; k < plot.size() / 3 * 2; k++) {
    plot[k] = 1;
  }
  for (Index k = plot.size() / 3 * 2; k < plot.size(); k++) {
    plot[k] = 0;
  }
  fft(plot.data(), plot_fft.data(), plot.size());
  NumpySaver("build/output/plot.npy")
      << plot << plot_fft.real() << plot_fft.imag();

  compress_image("images/A.png", "build/output/A9.ldw", 0.9);
  compress_image("images/A.png", "build/output/A5.ldw", 0.5);
  compress_image("images/A.png", "build/output/A2.ldw", 0.2);
  compress_image("images/A.png", "build/output/A1.ldw", 0.1);
  decompress_image("build/output/A9.ldw", "build/plots/A9.png");
  decompress_image("build/output/A5.ldw", "build/plots/A5.png");
  decompress_image("build/output/A2.ldw", "build/plots/A2.png");
  decompress_image("build/output/A1.ldw", "build/plots/A1.png");

  blur("images/dune.png", "build/output/dune_blur1.png",
       "build/output/dune_blur1_mask.png", .1);
  blur("images/dune.png", "build/output/dune_blur2.png",
       "build/output/dune_blur2_mask.png", .2);
  blur("images/dune.png", "build/output/dune_blur4.png",
       "build/output/dune_blur4_mask.png", .4);
  blur("images/dune.png", "build/output/dune_blur3.png",
       "build/output/dune_blur3_mask.png", .3);

  blur_smooth("images/dune.png", "build/output/dune_blur_smooth.png",
              "build/output/dune_blur_smooth_mask.png", .001);

  rect_filter("images/dune.png", "build/output/dune_blur_rect.png",
              "build/output/dune_rect_mask.png", .1, .1);

  sharpen("images/blackhole.png", "build/output/blackhole_sharp9.png",
          "build/output/blackhole_sharp9_mask.png", .09);
  sharpen("images/blackhole.png", "build/output/blackhole_sharp6.png",
          "build/output/blackhole_sharp6_mask.png", .06);
  sharpen("images/blackhole.png", "build/output/blackhole_sharp3.png",
          "build/output/blackhole_sharp3_mask.png", .03);
  sharpen("images/blackhole.png", "build/output/blackhole_sharp1.png",
          "build/output/blackhole_sharp1_mask.png", .01);

  sharpen_smooth("images/blackhole.png",
                 "build/output/blackhole_sharp_smooth_g.png",
                 "build/output/blackhole_sharp_smooth_mask_g.png", .01, true);
  sharpen_smooth("images/blackhole.png",
                 "build/output/blackhole_sharp_smooth.png",
                 "build/output/blackhole_sharp_smooth_mask.png", .01, false);

  // compress_image("images/dune.png", "build/output/dune8.ldw", 0.9);
  // compress_image("images/dune.png", "build/output/dune5.ldw", 0.5);
  // compress_image("images/dune.png", "build/output/dune1.ldw", 0.1);
  // decompress_image("build/output/dune5.ldw", "build/plots/dune9.png");
  // decompress_image("build/output/dune5.ldw", "build/plots/dune5.png");
  // decompress_image("build/output/dune1.ldw", "build/plots/dune1.png");

  //   sharpen("images/blackhole.png", "build/plots/blackhole001_bw.png",
  //           "build/plots/blackhole001_bw_mask.png", .001, true);
  //   sharpen("images/blackhole.png", "build/plots/blackhole05_bw.png",
  //           "build/plots/blackhole01_bw_mask.png", .01, true);
  //   sharpen("images/blackhole.png", "build/plots/blackhole01_bw.png",
  //           "build/plots/blackhole01_bw_mask.png", .01, true);
  //   sharpen("images/blackhole.png", "build/plots/blackhole0_bw.png",
  //           "build/plots/blackhole0_bw_mask.png", .0, true);

  //   sharpen("images/dune.png", "build/plots/dune_bw_sharp.png",
  //           "build/plots/dune_bw_sharp_mask.png", .005, true);
  //   blur("images/dune.png", "build/plots/dune_bw_blur.png",
  //        "build/plots/dune_bw_blur_mask.png", .95, true);
  //   sharpen("images/dune.png", "build/plots/dune_bw.png",
  //           "build/plots/dune_bw_mask.png", 0, true);

  //   rect_filter("images/dune.png", "build/plots/dune_bw_rect.png",
  //               "build/plots/dune_bw_rect_mask.png", .5, .5, true);
  //   anti_rect_filter("images/dune.png", "build/plots/dune_bw_arect.png",
  //                    "build/plots/dune_bw_arect_mask.png", .05, .05, true);

  //   blur("images/A.png", "build/plots/A_blur95.png",
  //        "build/plots/A_blur95_mask.png", .95);
  //   blur("images/A.png", "build/plots/A_blur9.png",
  //        "build/plots/A_blur9_mask.png", .9);
  //   blur("images/A.png", "build/plots/A_blur8.png",
  //        "build/plots/A_blur8_mask.png", .8);
  //   blur("images/A.png", "build/plots/A_blur7.png",
  //        "build/plots/A_blur7_mask.png", .7);
  //   blur("images/A.png", "build/plots/A_blur6.png",
  //        "build/plots/A_blur6_mask.png", .6);
  //   blur("images/A.png", "build/plots/A_blur5.png",
  //        "build/plots/A_blur5_mask.png", .5);
  //   blur("images/A.png", "build/plots/A_blur4.png",
  //        "build/plots/A_blur4_mask.png", .4);
  //   blur("images/A.png", "build/plots/A_blur3.png",
  //        "build/plots/A_blur3_mask.png", .3);
  //   blur("images/A.png", "build/plots/A_blur2.png",
  //        "build/plots/A_blur2_mask.png", .2);
  //   blur("images/A.png", "build/plots/A_blur1.png",
  //        "build/plots/A_blur1_mask.png", .1);

  //   sharpen("images/A.png", "build/plots/A_sharpen95.png",
  //           "build/plots/A_sharpen95_mask.png", .95);
  //   sharpen("images/A.png", "build/plots/A_sharpen9.png",
  //           "build/plots/A_sharpen9_mask.png", .9);
  //   sharpen("images/A.png", "build/plots/A_sharpen8.png",
  //           "build/plots/A_sharpen8_mask.png", .8);
  //   sharpen("images/A.png", "build/plots/A_sharpen7.png",
  //           "build/plots/A_sharpen7_mask.png", .7);
  //   sharpen("images/A.png", "build/plots/A_sharpen6.png",
  //           "build/plots/A_sharpen6_mask.png", .6);
  //   sharpen("images/A.png", "build/plots/A_sharpen5.png",
  //           "build/plots/A_sharpen5_mask.png", .5);
  //   sharpen("images/A.png", "build/plots/A_sharpen4.png",
  //           "build/plots/A_sharpen4_mask.png", .4);
  //   sharpen("images/A.png", "build/plots/A_sharpen3.png",
  //           "build/plots/A_sharpen3_mask.png", .3);
  //   sharpen("images/A.png", "build/plots/A_sharpen2.png",
  //           "build/plots/A_sharpen2_mask.png", .2);
  //   sharpen("images/A.png", "build/plots/A_sharpen1.png",
  //           "build/plots/A_sharpen1_mask.png", .1);

  //   rect_filter("images/A.png", "build/plots/A_blur_rect95.png",
  //               "build/plots/A_blur_rect95_mask.png", .95, .95);
  //   rect_filter("images/A.png", "build/plots/A_blur_rect9.png",
  //               "build/plots/A_blur_rect9_mask.png", .9, .9);
  //   rect_filter("images/A.png", "build/plots/A_blur_rect8.png",
  //               "build/plots/A_blur_rect8_mask.png", .8, .8);
  //   rect_filter("images/A.png", "build/plots/A_blur_rect7.png",
  //               "build/plots/A_blur_rect7_mask.png", .7, .7);
  //   rect_filter("images/A.png", "build/plots/A_blur_rect6.png",
  //               "build/plots/A_blur_rect6_mask.png", .6, .6);
  //   rect_filter("images/A.png", "build/plots/A_blur_rect5.png",
  //               "build/plots/A_blur_rect5_mask.png", .5, .5);
  //   rect_filter("images/A.png", "build/plots/A_blur_rect4.png",
  //               "build/plots/A_blur_rect4_mask.png", .4, .4);
  //   rect_filter("images/A.png", "build/plots/A_blur_rect3.png",
  //               "build/plots/A_blur_rect3_mask.png", .3, .3);
  //   rect_filter("images/A.png", "build/plots/A_blur_rect2.png",
  //               "build/plots/A_blur_rect2_mask.png", .2, .2);
  //   rect_filter("images/A.png", "build/plots/A_blur_rect1.png",
  //               "build/plots/A_blur_rect1_mask.png", .1, .1);

  //   anti_rect_filter("images/A.png", "build/plots/A_blur_anti_rect95.png",
  //                    "build/plots/A_blur_anti_rect95_mask.png", .95, .95);
  //   anti_rect_filter("images/A.png", "build/plots/A_blur_anti_rect9.png",
  //                    "build/plots/A_blur_anti_rect9_mask.png", .9, .9);
  //   anti_rect_filter("images/A.png", "build/plots/A_blur_anti_rect8.png",
  //                    "build/plots/A_blur_anti_rect8_mask.png", .8, .8);
  //   anti_rect_filter("images/A.png", "build/plots/A_blur_anti_rect7.png",
  //                    "build/plots/A_blur_anti_rect7_mask.png", .7, .7);
  //   anti_rect_filter("images/A.png", "build/plots/A_blur_anti_rect6.png",
  //                    "build/plots/A_blur_anti_rect6_mask.png", .6, .6);
  //   anti_rect_filter("images/A.png", "build/plots/A_blur_anti_rect5.png",
  //                    "build/plots/A_blur_anti_rect5_mask.png", .5, .5);
  //   anti_rect_filter("images/A.png", "build/plots/A_blur_anti_rect4.png",
  //                    "build/plots/A_blur_anti_rect4_mask.png", .4, .4);
  //   anti_rect_filter("images/A.png", "build/plots/A_blur_anti_rect3.png",
  //                    "build/plots/A_blur_anti_rect3_mask.png", .3, .3);
  //   anti_rect_filter("images/A.png", "build/plots/A_blur_anti_rect2.png",
  //                    "build/plots/A_blur_anti_rect2_mask.png", .2, .2);
  //   anti_rect_filter("images/A.png", "build/plots/A_blur_anti_rect1.png",
  //                    "build/plots/A_blur_anti_rect1_mask.png", .1, .1);

  //   sharpen("images/boy.png", "build/plots/boy.png",
  //           "build/plots/boy_control.png", .05);

  //   sharpen_smooth("images/boy.png", "build/plots/boy_smooth.png",
  //                  "build/plots/boy_control_smooth.png", .001);
}
