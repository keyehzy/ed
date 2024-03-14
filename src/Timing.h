#pragma once

#include <chrono>
#include <functional>

class TimeIt {
 public:
  TimeIt(std::function<void()> body);

  std::chrono::duration<double> duration() const { return m_duration; }

 private:
  std::chrono::duration<double> m_duration;
};
