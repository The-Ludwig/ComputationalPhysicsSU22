

#include <yaml-cpp/yaml.h>

#include <Eigen/Dense>
#include <iostream>

#include "NumpySaver.hpp"

int main(int argc, char const *argv[]) {
  YAML::Node config = YAML::LoadFile("config.yaml");

  Eigen::ArrayXd array;
  std::string name;

  if (config["dim"] && config["name"]) {
    array = Eigen::ArrayXd::LinSpaced(config["dim"].as<int>(), 0, 10);
    name = config["name"].as<std::string>();
  } else {
    throw std::range_error("Config file must include 'dim_x' and 'dim_y'.");
  }

  NumpySaver("build/output/" + name + ".npy")
      << "This is a test array" << array << (array * array * 2 + 10.);

  std::cout << "build/output/test.npy saved!";

  return 0;
}
