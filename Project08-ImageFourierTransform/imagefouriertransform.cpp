#include "image_fft.hpp"

int main(int argc, char const* argv[]) {
#ifndef NDEBUG
  FFT::tests();
#endif
  image_test();
}
