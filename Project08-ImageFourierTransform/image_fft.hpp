#pragma once

#include <Eigen/Dense>
#include <boost/gil.hpp>
#include <boost/gil/extension/io/png.hpp>
#include <cassert>
#include <iostream>
#include <vector>

#include "fft.hpp"
#include "progressbar.hpp"

using namespace boost::gil;
using Eigen::MatrixXcd;

template <typename GrayView>
Eigen::MatrixXcd dft_gray(const GrayView& src) {
  typedef std::complex<double> complex;
  typedef Eigen::Matrix<complex, -1, -1, Eigen::ColMajor> mattype;

  auto h = src.height();
  auto w = src.width();
  mattype dft(h, w);

  progressbar bar(h + w);

  // first do fourier traffo of rows
  for (int y = 0; y < h; y++) {
    bar.update();
    typename GrayView::x_iterator src_it = src.row_begin(y);

    std::vector<complex> buf;
    buf.resize(w);
    fft(src_it, buf.data(), w);

    for (int x = 0; x < w; x++) dft(y, x) = buf[x];
  }

  // now of cols
  for (int x = 0; x < w; x++) {
    bar.update();
    std::vector<complex> buf;
    buf.resize(h);
    fft(&dft(0, x), buf.data(), h);

    for (int y = 0; y < h; y++) dft(y, x) = buf[y];
  }

  return dft;
}

template <typename SrcGrayView, typename DstGrayView>
void naive_dft_gray(const SrcGrayView& src, DstGrayView& dst) {
  assert(src.dimensions() == dst.dimensions());

  MatrixXcd mat = dft_gray(src);

  auto w = dst.width();
  auto h = dst.height();

  for (int y = 0; y < dst.height(); y++) {
    typename DstGrayView::x_iterator dst_it = dst.row_begin(y);
    for (int x = 0; x < dst.width(); x++)
      dst_it[x] = abs(mat((y + h / 2) % h, (x + w / 2) % w));
  }
}

void image_test() {
  using namespace boost::gil;

  std::string filename("images/webb.png");
  rgb8_image_t img;
  read_image(filename, img, png_tag());
  std::cout << "Image read! Dimensions x: " << img.dimensions().x
            << " y: " << img.dimensions().y << std::endl;

  gray8_image_t img_ft(img.dimensions());
  std::cout << "Now performing a gray naive fourier transform" << std::endl;
  naive_dft_gray(color_converted_view<gray8_pixel_t>(view(img)), view(img_ft));

  std::cout << "Saving image" << std::endl;
  write_view("build/output/test.png", view(img_ft), png_tag());
}
