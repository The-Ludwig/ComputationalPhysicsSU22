#pragma once

#include <Eigen/Dense>
#include <boost/gil.hpp>
#include <boost/gil/extension/io/png.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <cassert>
#include <cinttypes>
#include <fstream>
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
  std::vector<double> buf;
  buf.resize(w);
  for (int y = 0; y < h; y++) {
    typename ImgView::x_iterator src_it = src.row_begin(y);
    for (size_t c = 0; c < nc; c++) {
      bar.update();
      auto src_channel_it =
          channel_iterator<cs_t, typename ImgView::x_iterator>(src_it, c);
      dct(src_channel_it, buf.data(), w);
      for (int x = 0; x < w; x++) dfts[c](y, x) = buf[x];
    }
  }

  // now of cols
  buf.resize(h);
  for (int x = 0; x < w; x++) {
    for (size_t c = 0; c < nc; c++) {
      bar.update();
      dct(&(dfts[c](0, x)), buf.data(), h);

      for (int y = 0; y < h; y++) dfts[c](y, x) = buf[y];
    }
  }

  return dfts;
}

template <typename ImgView>
vector<Eigen::MatrixXd> idct(const ImgView& src) {
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
      idct(src_channel_it, buf.data(), w);
      for (int x = 0; x < w; x++) dfts[c](y, x) = buf[x];
    }
  }

  // now of cols
  for (int x = 0; x < w; x++) {
    for (size_t c = 0; c < nc; c++) {
      bar.update();
      std::vector<double> buf;
      buf.resize(h);
      idct(&(dfts[c](0, x)), buf.data(), h);

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
  double max = std::numeric_limits<double>::min();
  double min = std::numeric_limits<double>::max();
  for (Eigen::Index row = 0; row < mat.rows(); row++)
    for (Eigen::Index col = 0; col < mat.cols(); col++) {
      double val = mat(row, col);
      if (val > max) max = val;
      if (val < min) min = val;
    }
  return {min, max};
}

std::tuple<double, double> real_bounds(MatrixXcd& mat) {
  double max = std::numeric_limits<double>::min();
  double min = std::numeric_limits<double>::max();
  for (Eigen::Index row = 0; row < mat.rows(); row++)
    for (Eigen::Index col = 0; col < mat.cols(); col++) {
      double val = mat(row, col).real();
      if (val > max) max = val;
      if (val < min) min = val;
    }
  return {min, max};
}

template <typename SrcView, typename DstView>
void dft(const SrcView& src, DstView& dst_mag, DstView& dst_phase,
         bool shifted = true, bool log = true) {
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
        Eigen::Index xcord, ycord;
        if (shifted) {
          ycord = (y + h / 2) % h;
          xcord = (x + w / 2) % w;
        } else {
          xcord = x;
          ycord = y;
        }
        if (log) {
          dst_mag_it[x][c] =
              (std::log(std::abs(mats[c](xcord, ycord))) - std::log(min)) /
              (std::log(max) - std::log(min)) * max_val;
        } else {
          dst_mag_it[x][c] = ((std::abs(mats[c](xcord, ycord))) - (min)) /
                             ((max) - (min)) * max_val;
        }
        dst_phase_it[x][c] =
            std::arg(mats[c](xcord, ycord)) / 2. / pi * max_val;
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
        dst_it[x][c] = mats[c](y, x);
      }
  }
}

template <typename SrcView, typename DstView>
void idct(const SrcView& src, DstView& dst) {
  assert(src.dimensions() == dst.dimensions());

  typedef typename channel_type<DstView>::type cs_t;
  cs_t max_val = std::numeric_limits<cs_t>::max();

  auto mats = idct(src);

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

size_t closest_smaller_power2(size_t number) {
  size_t log2n = 0;
  while ((number >> ++log2n) > 0) {
  };

  return 1 << (log2n - 1);
}

template <typename DstView>
void to_image(vector<MatrixXcd>& src, DstView& dst) {
  assert(src[0].cols() == dst.width());
  assert(src[0].rows() == dst.height());

  typedef typename channel_type<DstView>::type cs_t;
  cs_t max_val = std::numeric_limits<cs_t>::max();

  std::vector<std::tuple<double, double>> limits;
  limits.reserve(src.size());
  for (size_t k = 0; k < src.size(); k++) limits.push_back(real_bounds(src[k]));

  auto w = src[0].cols();
  auto h = src[0].rows();

  for (int y = 0; y < h; y++) {
    typename DstView::x_iterator dst_it = dst.row_begin(y);
    for (int x = 0; x < w; x++)
      for (size_t c = 0; c < src.size(); c++) {
        auto [min, max] = limits[c];

        dst_it[x][c] = src[c](y, x).real();
      }
  }
}

void apply_filter(const char* from, const char* to,
                  std::function<void(MatrixXcd&)> filter) {
  rgb8_image_t img;
  read_image(from, img, png_tag());

  typedef typename channel_type<rgb8_image_t>::type cs_t;
  cs_t max_val = std::numeric_limits<cs_t>::max();

  auto mats = dft(view(img));

  auto w = img.width();
  auto h = img.height();

  // apply filter
  for (size_t c = 0; c < mats.size(); c++) filter(mats[c]);

  // transform back
  // first do fourier traffo of rows
  progressbar bar((h + w) * mats.size());
  std::vector<complex<double>> buf;
  buf.resize(w);
  for (int y = 0; y < h; y++) {
    for (size_t c = 0; c < mats.size(); c++) {
      bar.update();
      ifft(mats[c].row(y), buf.data(), w);
      for (int x = 0; x < w; x++) mats[c](y, x) = buf[x];
    }
  }
  // now of columns
  buf.resize(h);
  for (int x = 0; x < w; x++) {
    for (size_t c = 0; c < mats.size(); c++) {
      bar.update();
      ifft(mats[c].col(x), buf.data(), h);
      for (int y = 0; y < h; y++) mats[c](y, x) = buf[y];
    }
  }

  rgb8_image_t sharpened(img.dimensions());
  to_image(mats, view(sharpened));

  write_view(to, view(sharpened), png_tag());
}

void sharpen(const char* from, const char* to, double radius) {
  apply_filter(from, to, [&](MatrixXcd& mat) {
    double limit = radius * std::min(mat.cols(), mat.rows()) / 2.;
    limit *= limit;
    for (int x = 0; x < mat.cols(); x++)
      for (int y = 0; y < mat.rows(); y++) {
        double pos_x =
            ((x + mat.cols() / 2) % mat.cols()) - (mat.cols() - 1) / 2.;
        double pos_y =
            ((y + mat.rows() / 2) % mat.rows()) - (mat.rows() - 1) / 2.;
        if (pos_x * pos_x + pos_y * pos_y < limit) {
          mat(y, x) = 0;
        }
      }
  });
}

void blur(const char* from, const char* to, double radius) {
  apply_filter(from, to, [&](MatrixXcd& mat) {
    double limit = radius * std::min(mat.cols(), mat.rows()) / 2.;
    for (int x = 0; x < mat.cols(); x++)
      for (int y = 0; y < mat.rows(); y++) {
        double pos_x =
            ((x + mat.cols() / 2) % mat.cols()) - (mat.cols() - 1) / 2.;
        double pos_y =
            ((y + mat.rows() / 2) % mat.rows()) - (mat.rows() - 1) / 2.;
        if (pos_x * pos_x + pos_y * pos_y > limit * limit) {
          mat(y, x) = 0;
        }
      }
  });
}

typedef uint_least16_t pos_t;
typedef uint_least8_t channel_t;
typedef std::tuple<pos_t, pos_t, channel_t> index_t;
void compress_image(const char* from, const char* to,
                    double compression_level) {
  rgb8_image_t img;
  read_image(from, img, png_tag());
  std::cout << "Image read! Dimensions x: " << img.dimensions().x
            << " y: " << img.dimensions().y << std::endl;

  typedef typename channel_type<rgb8_image_t>::type cs_t;
  cs_t max_val = std::numeric_limits<cs_t>::max();

  auto mats = dct(view(img));
  // x, y, c
  std::vector<index_t> indices;

  indices.reserve(mats[0].cols() * mats[0].rows() * mats.size());
  for (size_t x = 0; x < mats[0].cols(); x++)
    for (size_t y = 0; y < mats[0].rows(); y++)
      for (size_t c = 0; c < mats.size(); c++) indices.push_back({x, y, c});

  std::sort(indices.begin(), indices.end(),
            [&](index_t& a, index_t& b) -> bool {
              auto [xa, ya, ca] = a;
              auto [xb, yb, cb] = b;
              return std::abs(mats[ca](ya, xa)) > std::abs(mats[cb](yb, xb));
            });

  std::fstream file(to, file.binary | file.trunc | file.out);
  if (!file.is_open()) {
    throw std::runtime_error("failed to open file");
  } else {
    boost::iostreams::filtering_ostream os;
    os.push(boost::iostreams::zlib_compressor());
    os.push(file);
    // write dimensions
    uint_least16_t width = mats[0].cols();
    uint_least16_t height = mats[0].rows();
    os.write(reinterpret_cast<char*>(&width), sizeof width);  // binary output
    os.write(reinterpret_cast<char*>(&height),
             sizeof height);  // binary output

    size_t lastk = indices.size() * compression_level;
    progressbar b2(lastk);

    for (size_t k = 0; k < lastk; k++) {
      b2.update();
      auto [x, y, c] = indices[k];
      float val = mats[c](y, x);

      os.write(reinterpret_cast<char*>(&x), sizeof x);      // binary output
      os.write(reinterpret_cast<char*>(&y), sizeof y);      // binary output
      os.write(reinterpret_cast<char*>(&c), sizeof c);      // binary output
      os.write(reinterpret_cast<char*>(&val), sizeof val);  // binary output
    }
  }
}

void decompress_image(const char* from, const char* to) {
  std::fstream file(from, file.binary | file.in);
  if (!file.is_open()) {
    throw std::runtime_error("failed to open file");
  }
  boost::iostreams::filtering_istream is;
  is.push(boost::iostreams::zlib_decompressor());
  is.push(file);
  // read dimensions
  uint_least16_t width, height;
  is.read(reinterpret_cast<char*>(&width), sizeof width);
  is.read(reinterpret_cast<char*>(&height), sizeof height);

  std::vector<MatrixXd> mats;
  for (int c = 0; c < 3; c++) mats.push_back(MatrixXd::Zero(height, width));

  while (!is.eof()) {
    pos_t x, y;
    channel_t c;
    float val;

    is.read(reinterpret_cast<char*>(&x), sizeof x);      // binary input
    is.read(reinterpret_cast<char*>(&y), sizeof y);      // binary input
    is.read(reinterpret_cast<char*>(&c), sizeof c);      // binary input
    is.read(reinterpret_cast<char*>(&val), sizeof val);  // binary input

    mats[c](y, x) = val;
  }

  // first do fourier traffo of rows
  std::vector<double> buf;
  buf.resize(width);
  for (int y = 0; y < height; y++) {
    for (size_t c = 0; c < 3; c++) {
      idct(mats[c].row(y), buf.data(), width);
      for (int x = 0; x < width; x++) mats[c](y, x) = buf[x];
    }
  }

  rgb8_image_t img(width, height);
  typedef typename channel_type<rgb8_image_t>::type cs_t;
  cs_t max_val = std::numeric_limits<cs_t>::max();
  auto img_view = view(img);

  // now of cols
  buf.resize(height);
  for (int x = 0; x < width; x++) {
    for (size_t c = 0; c < 3; c++) {
      idct(mats[c].col(x), buf.data(), height);
      for (int y = 0; y < height; y++) img_view(x, y)[c] = buf[y];
    }
  }

  write_view(to, img_view, png_tag());
}

void dft_image(const char* from, const char* to_mag, const char* to_phase,
               bool crop = true, bool log = true, bool shifted = true) {
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
    dft(sub_img, view(img_ft_mag), view(img_ft_phase), shifted, log);
  } else {
    img_ft_mag = rgb8_image_t(img.dimensions());
    img_ft_phase = rgb8_image_t(img.dimensions());
    dft(view(img), view(img_ft_mag), view(img_ft_phase), shifted, log);
  }

  write_view(to_mag, view(img_ft_mag), png_tag());
  write_view(to_phase, view(img_ft_phase), png_tag());
}

void dct_image(const char* from, const char* to, bool crop = false) {
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

void idct_image(const char* from, const char* to, bool crop = false) {
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
    idct(sub_img, view(img_dct));
  } else {
    img_dct = rgb8_image_t(img.dimensions());
    idct(view(img), view(img_dct));
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

  // dft(color_converted_view<gray8_pixel_t>(view(img)), view(img_ft_mag),
  // view(img_ft_phase));
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

void sharpen_test() {
  sharpen("images/dune.png", "build/output/dune_sharp.png", .1);
  blur("images/dune.png", "build/output/dune_blur.png", .5);
}

void compress_test() {
  compress_image("images/dune.png", "build/output/dune_cmp.ldw", .01);
  decompress_image("build/output/dune_cmp.ldw", "build/output/reconst.png");
}