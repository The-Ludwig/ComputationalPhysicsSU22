#pragma once

#include <chrono>

class Timer {
 private:
  std::chrono::steady_clock::time_point start;
  std::chrono::duration<double> elapsed_seconds;

 public:
  Timer() { reset(); };

  void reset() { start = std::chrono::steady_clock::now(); }

  double elapsed_s() {
    elapsed_seconds = std::chrono::steady_clock::now() - start;
    return elapsed_seconds.count();
  }

  double elapsed_ms() { return elapsed_s() * 1e3; }

  /**
   * @brief Tests how long a function executes on average in ms.
   *
   * Using the Welfords online algorithm to calculate the mean and standard
   * deviation.
   *
   * @param precision Time limit for measuring performance in ms
   * @return std::tuple<double, double> Mean and Standard deviation in ms
   */
  static std::tuple<double, double> measure_time(std::function<void()> function,
                                                 double precision = .1) {
    unsigned int n = 0;
    double mean = 0;
    double S = 0;

    Timer lap;
    do {
      lap.reset();
      function();
      double lap_time = lap.elapsed_ms();
      n++;
      double old_mean = mean;
      mean += (lap_time - mean) / double(n);
      S += (lap_time - mean) * (lap_time - old_mean);
    } while (n < 10 || std::sqrt(S / ((n - 1) * n)) > precision);

    return {mean, std::sqrt(S / (n - 1) / n)};
  }
};