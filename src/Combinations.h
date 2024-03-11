#pragma once

#include "Assert.h"

// Generate all possible combinations of states with size S and filling N
// Ex: for S = 2, N = 2
// There are 6 states with size 4 and 2 bits set
// 0011 0101 1001 0110 1010 1100

template <size_t S, size_t N>
class Combinations {
 public:
  Combinations() { generate_combinations(); }

  const std::vector<std::array<uint8_t, S>>& states() const { return m_states; }

  size_t size() const { return m_states.size(); }

  std::array<uint8_t, S> i_to_b(size_t i) const {
    ED_ASSERT(!m_states.empty());
    ED_ASSERT(i < m_states.size());
    return m_states[i];
  }

  int b_to_i(const std::array<uint8_t, S>& b) const {
    ED_ASSERT(!m_states.empty());
    auto it = std::find(m_states.begin(), m_states.end(), b);
    ED_ASSERT(it != m_states.end());
    return std::distance(m_states.begin(), it);
  }

 private:
  void generate_combinations() {
    if (!m_states.empty()) {
      return;
    }
    std::array<uint8_t, S> v{};
    std::fill(v.begin() + S - N, v.end(), 1);
    do {
      m_states.emplace_back(v);
    } while (std::next_permutation(v.begin(), v.end()));
  }

  std::vector<std::array<uint8_t, S>> m_states;
};
