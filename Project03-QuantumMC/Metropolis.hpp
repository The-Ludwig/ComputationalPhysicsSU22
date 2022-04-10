#pragma once

#include <Eigen/Dense>
#include <random>

template <int NumberArguments>
class MetropolisAlgorithm {
  typedef Array<double, NumArguments, 1> Argument;
  typedef double(const Argument&) FunctionType;
  typedef std::mt19937 Generator;

 private:
  std::function<FunctionType> pdf;
  std::function<double(Generator)> step;
  std::function<double(Generator)> uniform;
  std::mt19937 gen;

  Argument current;

 public:
  MetropolisAlgorithm(std::function<pdf> pdf, Argument argstart,
                      std::function<double(Generator)> step,
                      std::function<double(Generator)> uniform,
                      Generator::return_type seed) {}
};