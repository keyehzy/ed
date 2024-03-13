#pragma once

#include <chrono>

class TimeIt {
 public:
  TimeIt(std::function<void()> body) {
    auto start = std::chrono::steady_clock::now();
    body();
    auto end = std::chrono::steady_clock::now();
    m_duration = end - start;
  }

  std::chrono::duration<double> duration() const { return m_duration; }

 private:
  std::chrono::duration<double> m_duration;
};
