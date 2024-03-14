#include "Timing.h"

TimeIt::TimeIt(std::function<void()> body) {
  auto start = std::chrono::steady_clock::now();
  body();
  auto end = std::chrono::steady_clock::now();
  m_duration = end - start;
}
