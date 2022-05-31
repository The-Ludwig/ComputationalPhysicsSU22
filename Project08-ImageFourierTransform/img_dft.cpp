#include <algorithm>
#include <iostream>

#include "image_fft.hpp"

int main(int argc, char const* argv[]) {
  std::vector<std::string> args(argv, argv + argc);

  if (argc < 4) {
    std::cout << "Please provide at least 3 arguments" << std::endl;
    return 0;
  }

  std::cout << "Opening '" << argv[1] << "'" << std::endl;

  bool log = std::find(args.begin(), args.end(), "--log") != args.end();
  bool crop = std::find(args.begin(), args.end(), "--crop") != args.end();
  bool shift = std::find(args.begin(), args.end(), "--noshift") == args.end();

  dft_image(argv[1], argv[2], argv[3], crop, log, shift);
}