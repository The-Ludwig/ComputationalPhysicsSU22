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
   * @brief Tests how long a function executes on average.
   *
   * Using the Welfords online algorithm to calculate the mean and standard
   * deviation.
   *
   * @param time_limit_s Time limit for measuring performance
   * @return std::tuple<double, double> Mean and Standard deviation in ms
   */
  static std::tuple<double, double> measure_time(std::function<void()> function,
                                                 double time_limit_s = 5) {
    unsigned int n = 0;
    double mean = 0;
    double mn = 0;

    Timer timer;
    Timer lap;
    while (timer.elapsed_s() < time_limit_s) {
      lap.reset();
      function();
      double lap_time = lap.elapsed_ms();
      n++;
      double add_mn = (lap_time - mean);
      mean += (lap_time - mean) / double(n);
      mn += (lap_time - mean) * add_mn;
    }

    return {mean, mn / (n - 1)};
  }
};