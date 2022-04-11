#pragma once

#include <Eigen/Dense>
#include <random>

template <int NumArguments, typename StepRNG = std::uniform_real_distribution<>,
          class Generator = std::mt19937>
class MetropolisAlgorithm {
 public:
  using Argument = Eigen::Array<double, NumArguments, 1>;
  using FunctionType = std::function<double(const Argument&)>;
  using GeneratorResType = Generator::result_type;

 private:
  FunctionType pdf;
  StepRNG step;
  std::uniform_real_distribution<> uniform =
      std::uniform_real_distribution<>(0, 1);
  std::mt19937 gen;

  double p_old;
  std::size_t argument_index = 0;

 public:
  Argument argument;

  MetropolisAlgorithm(FunctionType pdf, Argument& argstart, StepRNG step,
                      GeneratorResType seed)
      : pdf(pdf), argument(argstart), step(step) {
    gen = Generator(seed);  // Standard mersenne_twister_engine seeded with seed

    p_old = pdf(argument);
  };

  MetropolisAlgorithm(FunctionType pdf, Argument& argstart, StepRNG step)
      : pdf(pdf), argument(argstart), step(step) {
    std::random_device rd;
    gen = Generator(rd());  // Standard mersenne_twister_engine seeded with seed

    p_old = pdf(argument);
  };

  void do_step() {
    double x_old = argument[argument_index];
    double x_new = x_old + step(gen);

    argument[argument_index] = x_new;
    double p_new = pdf(argument);
    double p = p_new / p_old;

    if (p < 1) {
      double s = uniform(gen);

      // reject new values
      if (s > p) {
        x_new = x_old;
        p_new = p_old;
      }
    }

    argument[argument_index] = x_new;
    p_old = p_new;

    argument_index++;
    if (argument_index >= NumArguments) {
      argument_index = 0;
    }
  }

  Eigen::Array<double, Eigen::Dynamic, NumArguments> get_sample(
      Eigen::Index N) {
    using namespace Eigen;
    typedef Array<double, Eigen::Dynamic, NumArguments> ReturnType;

    ReturnType returner(N);

    for (Index i = 0; i < N; i++) {
      do_step();
      returner.row(i) = argument;
    }

    return returner;
  }
};