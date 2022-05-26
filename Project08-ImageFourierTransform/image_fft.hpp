#pragma once

#include <Eigen/Dense>
#include <boost/gil.hpp>
#include <boost/gil/extension/io/png.hpp>
#include <cassert>
#include <iostream>
#include <limits>  // std::numeric_limits
#include <typeinfo>
#include <vector>

#include "fft.hpp"
#include "progressbar.hpp"

using namespace boost::gil;
using Eigen::MatrixXcd, Eigen::MatrixXd;
using std::vector;

// member typedefs provided through inheriting from std::iterator
template <typename Pixel_t, typename Original_t>
class channel_iterator
    : public std::iterator<
          std::random_access_iterator_tag,  // iterator_category
          Pixel_t                           // value_type
          > {
  const Original_t& original;
  size_t channel;

 public:
  explicit channel_iterator(const Original_t& original, size_t channel)
      : original(original), channel(channel){};
  const Pixel_t& operator[](size_t n) const { return original[n][channel]; };
};

template <typename ImgView>
vector<Eigen::MatrixXcd> dft(const ImgView& src) {
  typedef std::complex<double> complex;
  typedef Eigen::Matrix<complex, -1, -1, Eigen::ColMajor> mattype;

  typedef typename channel_type<ImgView>::type cs_t;

  auto h = src.height();
  auto w = src.width();
  auto nc = num_channels<ImgView>::value;
  vector<mattype> dfts;
  dfts.reserve(nc);
  for (size_t k = 0; k < nc; k++) dfts.push_back(mattype(h, w));

  progressbar bar((h + w) * nc);

  // first do fourier traffo of rows
  for (int y = 0; y < h; y++) {
    typename ImgView::x_iterator src_it = src.row_begin(y);
    for (size_t c = 0; c < nc; c++) {
      bar.update();
      auto src_channel_it =
          channel_iterator<cs_t, typename ImgView::x_iterator>(src_it, c);
      std::vector<complex> buf;
      buf.resize(w);
      fft(src_channel_it, buf.data(), w);
      for (int x = 0; x < w; x++) dfts[c](y, x) = buf[x];
    }
  }

  // now of cols
  for (int x = 0; x < w; x++) {
    for (size_t c = 0; c < nc; c++) {
      bar.update();
      std::vector<complex> buf;
      buf.resize(h);
      fft(&(dfts[c](0, x)), buf.data(), h);

      for (int y = 0; y < h; y++) dfts[c](y, x) = buf[y];
    }
  }

  return dfts;
}

template <typename ImgView>
vector<Eigen::MatrixXd> dct(const ImgView& src) {
  typedef Eigen::Matrix<double, -1, -1, Eigen::ColMajor> mattype;

  typedef typename channel_type<ImgView>::type cs_t;

  auto h = src.height();
  auto w = src.width();
  auto nc = num_channels<ImgView>::value;
  vector<mattype> dfts;
  dfts.reserve(nc);
  for (size_t k = 0; k < nc; k++) dfts.push_back(mattype(h, w));

  progressbar bar((h + w) * nc);

  // first do fourier traffo of rows
  for (int y = 0; y < h; y++) {
    typename ImgView::x_iterator src_it = src.row_begin(y);
    for (size_t c = 0; c < nc; c++) {
      bar.update();
      auto src_channel_it =
          channel_iterator<cs_t, typename ImgView::x_iterator>(src_it, c);
      std::vector<double> buf;
      buf.resize(w);
      dct(src_channel_it, buf.data(), w);
      for (int x = 0; x < w; x++) dfts[c](y, x) = buf[x];
    }
  }

  // now of cols
  for (int x = 0; x < w; x++) {
    for (size_t c = 0; c < nc; c++) {
      bar.update();
      std::vector<double> buf;
      buf.resize(h);
      dct(&(dfts[c](0, x)), buf.data(), h);

      for (int y = 0; y < h; y++) dfts[c](y, x) = buf[y];
    }
  }

  return dfts;
}

std::tuple<double, double> mag_bounds(MatrixXcd& mat) {
  double max = 0;
  double min = std::numeric_limits<double>::max();
  for (Eigen::Index row = 0; row < mat.rows(); row++)
    for (Eigen::Index col = 0; col < mat.cols(); col++) {
      double val = std::abs(mat(row, col));
      if (val > max) max = val;
      if (val < min) min = val;
    }
  return {min, max};
}

std::tuple<double, double> bounds(MatrixXd& mat) {
  double max = 0;
  double min = std::numeric_limits<double>::max();
  for (Eigen::Index row = 0; row < mat.rows(); row++)
    for (Eigen::Index col = 0; col < mat.cols(); col++) {
      double val = mat(row, col);
      if (val > max) max = val;
      if (val < min) min = val;
    }
  return {min, max};
}

template <typename SrcView, typename DstView>
void dft(const SrcView& src, DstView& dst_mag, DstView& dst_phase,
         bool shifted = true) {
  assert(src.dimensions() == dst_mag.dimensions());
  assert(src.dimensions() == dst_phase.dimensions());

  typedef typename channel_type<DstView>::type cs_t;
  cs_t max_val = std::numeric_limits<cs_t>::max();

  auto mats = dft(src);

  std::vector<std::tuple<double, double>> limits;
  limits.reserve(mats.size());
  for (size_t k = 0; k < mats.size(); k++)
    limits.push_back(mag_bounds(mats[k]));

  auto w = dst_mag.width();
  auto h = dst_mag.height();

  for (int y = 0; y < h; y++) {
    typename DstView::x_iterator dst_mag_it = dst_mag.row_begin(y);
    typename DstView::x_iterator dst_phase_it = dst_phase.row_begin(y);
    for (int x = 0; x < w; x++)
      for (size_t c = 0; c < mats.size(); c++) {
        auto [min, max] = limits[c];
        if (shifted) {
          dst_mag_it[x][c] =
              (std::log(std::abs(mats[c]((y + h / 2) % h, (x + w / 2) % w))) -
               std::log(min)) /
              (std::log(max) - std::log(min)) * max_val;
          dst_phase_it[x][c] =
              std::arg(mats[c]((y + h / 2) % h, (x + w / 2) % w)) / 2. / pi *
              max_val;
        } else {
          dst_mag_it[x][c] =
              (std::abs(mats[c](y, x)) - min) / (max - min) * max_val;
          ;
          dst_phase_it[x][c] = std::arg(mats[c](y, x)) / 2. / pi * max_val;
        }
      }
  }
}

template <typename SrcView, typename DstView>
void dct(const SrcView& src, DstView& dst) {
  assert(src.dimensions() == dst.dimensions());

  typedef typename channel_type<DstView>::type cs_t;
  cs_t max_val = std::numeric_limits<cs_t>::max();

  auto mats = dct(src);

  std::vector<std::tuple<double, double>> limits;
  limits.reserve(mats.size());
  for (size_t k = 0; k < mats.size(); k++) limits.push_back(bounds(mats[k]));

  auto w = dst.width();
  auto h = dst.height();

  for (int y = 0; y < h; y++) {
    typename DstView::x_iterator dst_it = dst.row_begin(y);
    for (int x = 0; x < w; x++)
      for (size_t c = 0; c < mats.size(); c++) {
        auto [min, max] = limits[c];
        dst_it[x][c] = (mats[c](y, x) - min) / (max - min) * max_val;
        ;
      }
  }
}

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
void dft_gray(const SrcGrayView& src, DstGrayView& dst_mag,
              DstGrayView& dst_phase, bool shifted = true) {
  assert(src.dimensions() == dst_mag.dimensions());
  assert(src.dimensions() == dst_phase.dimensions());

  MatrixXcd mat = dft_gray(src);

  auto w = dst_mag.width();
  auto h = dst_mag.height();

  for (int y = 0; y < h; y++) {
    typename DstGrayView::x_iterator dst_mag_it = dst_mag.row_begin(y);
    typename DstGrayView::x_iterator dst_phase_it = dst_phase.row_begin(y);
    for (int x = 0; x < w; x++) {
      if (shifted) {
        dst_mag_it[x] = std::abs(mat((y + h / 2) % h, (x + w / 2) % w));
        dst_phase_it[x] = std::arg(mat((y + h / 2) % h, (x + w / 2) % w));
      } else {
        dst_mag_it[x] = std::abs(mat(y, x));
        dst_phase_it[x] = std::arg(mat(y, x));
      }
    }
  }
}

size_t closest_smaller_power2(size_t number) {
  size_t log2n = 0;
  while ((number >> ++log2n) > 0) {
  };

  return 1 << (log2n - 1);
}

void dft_image(const char* from, const char* to_mag, const char* to_phase,
               bool crop = true) {
  rgb8_image_t img;
  read_image(from, img, png_tag());
  std::cout << "Image read! Dimensions x: " << img.dimensions().x
            << " y: " << img.dimensions().y << std::endl;

  rgb8_image_t img_ft_mag, img_ft_phase;
  if (crop) {
    auto subw = closest_smaller_power2(img.width());
    auto subh = closest_smaller_power2(img.height());
    auto sub_img = subimage_view(view(img), 0, 0, subw, subh);
    img_ft_mag = rgb8_image_t(sub_img.dimensions());
    img_ft_phase = rgb8_image_t(sub_img.dimensions());
    dft(sub_img, view(img_ft_mag), view(img_ft_phase));
  } else {
    img_ft_mag = rgb8_image_t(img.dimensions());
    img_ft_phase = rgb8_image_t(img.dimensions());
    dft(view(img), view(img_ft_mag), view(img_ft_phase));
  }

  write_view(to_mag, view(img_ft_mag), png_tag());
  write_view(to_phase, view(img_ft_phase), png_tag());
}

void dct_image(const char* from, const char* to, bool crop = true) {
  rgb8_image_t img;
  read_image(from, img, png_tag());
  std::cout << "Image read! Dimensions x: " << img.dimensions().x
            << " y: " << img.dimensions().y << std::endl;

  rgb8_image_t img_dct;
  if (crop) {
    auto subw = closest_smaller_power2(img.width());
    auto subh = closest_smaller_power2(img.height());
    auto sub_img = subimage_view(view(img), 0, 0, subw, subh);
    img_dct = rgb8_image_t(sub_img.dimensions());
    dct(sub_img, view(img_dct));
  } else {
    img_dct = rgb8_image_t(img.dimensions());
    dct(view(img), view(img_dct));
  }

  write_view(to, view(img_dct), png_tag());
}

void image_test() {
  using namespace boost::gil;

  std::string filename("images/webb.png");
  rgb8_image_t img;
  read_image(filename, img, png_tag());
  std::cout << "Image read! Dimensions x: " << img.dimensions().x
            << " y: " << img.dimensions().y << std::endl;

  gray8_image_t img_ft_mag(img.dimensions());
  gray8_image_t img_ft_phase(img.dimensions());
  std::cout << "Now performing a gray naive fourier transform" << std::endl;

  dft(color_converted_view<gray8_pixel_t>(view(img)), view(img_ft_mag),
      view(img_ft_phase));
  std::cout << "Saving image" << std::endl;
  write_view("build/output/test_mag.png", view(img_ft_mag), png_tag());
  write_view("build/output/test_phase.png", view(img_ft_phase), png_tag());

  std::cout
      << "Now performing a fourier transform of closest power of two subimage"
      << std::endl;
  auto subw = closest_smaller_power2(img.width());
  auto subh = closest_smaller_power2(img.height());
  auto sub_img = subimage_view(color_converted_view<gray8_pixel_t>(view(img)),
                               0, 0, subw, subh);
  gray8_image_t img_sub_ft_mag(sub_img.dimensions());
  gray8_image_t img_sub_ft_phase(sub_img.dimensions());

  dft(sub_img, view(img_sub_ft_mag), view(img_sub_ft_phase));
  std::cout << "Saving image" << std::endl;
  write_view("build/output/sub_test_mag.png", view(img_sub_ft_mag), png_tag());
  write_view("build/output/sub_test_phase.png", view(img_sub_ft_phase),
             png_tag());

  std::cout << "Now performing a colored fourier transform of closest power of "
               "two subimage"
            << std::endl;
  auto sub_img_c = subimage_view(view(img), 0, 0, subw, subh);
  rgb8_image_t img_sub_c_ft_mag(sub_img_c.dimensions());
  rgb8_image_t img_sub_c_ft_phase(sub_img_c.dimensions());

  dft(sub_img_c, view(img_sub_c_ft_mag), view(img_sub_c_ft_phase));
  std::cout << "Saving image" << std::endl;
  write_view("build/output/sub_c_test_mag.png", view(img_sub_c_ft_mag),
             png_tag());
  write_view("build/output/sub_c_test_phase.png", view(img_sub_c_ft_phase),
             png_tag());

  std::cout << "Now performing a colored cosine transform of closest power of "
               "two subimage"
            << std::endl;
  rgb8_image_t img_sub_c_dct(sub_img_c.dimensions());

  dct(sub_img_c, view(img_sub_c_dct));
  std::cout << "Saving image" << std::endl;
  write_view("build/output/sub_c_test_dct.png", view(img_sub_c_dct), png_tag());
}
