#include <Eigen/Dense>
#include <boost/gil.hpp>
#include <boost/gil/extension/io/png.hpp>
#include <cassert>
#include <iostream>
#include <vector>

using std::vector;
using namespace boost::gil;
using Eigen::MatrixXcd;
using std::sqrt, std::exp, std::abs;
constexpr std::complex<double> i(0, 1);
constexpr double pi = 3.14159265358979323846;

template <typename GrayView>
Eigen::MatrixXcd naive_dft_gray(const GrayView& src) {
  typedef std::complex<double> complex;

  auto h = src.height();
  auto w = src.width();
  MatrixXcd dft(h, w);

  for (int y_ = 0; y_ < h; ++y_) {
    std::cout << "Row " << y_ << std::endl;
    for (int x_ = 0; x_ < w; ++x_) {
      complex f = 0;

      for (int y = 0; y < h; y++) {
        typename GrayView::x_iterator src_it = src.row_begin(y);
        for (int x = 0; x < w; x++)
          f += double(src_it[x][0]) *
               exp(-i * 2. * pi * (x * x_ / double(w) + y * y_ / double(h)));
      }

      dft(y_, x_) = f / sqrt(w * h);
    }
  }

  return dft;
}

template <typename SrcGrayView, typename DstGrayView>
void naive_dft_gray(const SrcGrayView& src, DstGrayView& dst) {
  MatrixXcd mat = naive_dft_gray(src);

  assert(src.dimensions() == dst.dimensions());

  for (int y = 0; y < dst.height(); y++) {
    typename DstGrayView::x_iterator dst_it = dst.row_begin(y);
    for (int x = 0; x < dst.width(); x++) dst_it[x] = abs(mat(y, x));
  }
}

int main(int argc, char const* argv[]) {
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
