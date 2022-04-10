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
};