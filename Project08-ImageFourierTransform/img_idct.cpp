#include <iostream>

#include "image_fft.hpp"

int main(int argc, char const* argv[]) {
  std::vector<std::string> args(argv, argv + argc);

  if (argc != 3 && argc != 4) {
    std::cout << "Please provide 2 or 3 arguments" << std::endl;
    return 0;
  }

  std::cout << "Opening '" << argv[1] << "'" << std::endl;

  if (argc == 3)
    idct_image(argv[1], argv[2]);
  else {
    if (args[3] == "false") idct_image(argv[1], argv[2], false);
    if (args[3] == "true") idct_image(argv[1], argv[2], true);
  }
}