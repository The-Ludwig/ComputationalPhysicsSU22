#include <iostream>

#include "image_fft.hpp"

int main(int argc, char const* argv[]) {
  std::vector<std::string> args(argv, argv + argc);

  if (argc != 4 && argc != 5) {
    std::cout << "Please provide 3 or 4 arguments" << std::endl;
    return 0;
  }

  std::cout << "Opening '" << argv[1] << "'" << std::endl;

  if (argc == 4)
    dft_image(argv[1], argv[2], argv[3]);
  else {
    if (args[4] == "false") dft_image(argv[1], argv[2], argv[3], false);
    if (args[4] == "true") dft_image(argv[1], argv[2], argv[3], true);
  }
}