#include <iostream>
#include <Eigen/Dense>
#include <yaml-cpp/yaml.h>
#include "NumpySaver.hpp"

int main(int argc, char const *argv[])
{

    YAML::Node config = YAML::LoadFile("config.yaml");

    Eigen::MatrixXd matrix;

    if (config["dim_x"] && config["dim_y"])
    {
        matrix = Eigen::MatrixXd::Ones(config["dim_x"].as<int>(), config["dim_y"].as<int>());
    }
    else
    {
        throw std::range_error("Config file must include 'dim_x' and 'dim_y'.");
    }

    std::cout << "Hello World!";

    NumpySaver("build/test.npy") << "Hello again World!" << matrix << matrix;

    std::cout << "build/test.npy saved!";

    return 0;
}
