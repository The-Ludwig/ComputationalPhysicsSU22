#include <boost/gil.hpp>
#include <boost/gil/extension/io/png.hpp>
#include <iostream>

template <typename GrayView>
void fourier_transform_gray(GrayView& img_src, GrayView& img_dest) {}

int main(int argc, char const* argv[]) {
  using namespace boost::gil;
  std::string filename("images/webb.png");
  rgb8_image_t img;
  read_image(filename, img, png_tag());

  std::cout << "x: " << img.dimensions().x << "y: " << img.dimensions().y;
}
