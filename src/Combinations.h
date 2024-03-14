#pragma once

#include <cstdint>
#include <vector>

#include "Assert.h"

// Generate all possible combinations of states with size S and filling N
// Ex: for S = 2, N = 2
// There are 6 states with size 4 and 2 bits set
// 0011 0101 1001 0110 1010 1100

template <size_t S, size_t N>
class Combinations {
 public:
  Combinations();

  const std::vector<std::array<uint8_t, S>>& states() const;

  size_t size() const;

  std::array<uint8_t, S> i_to_b(size_t i) const;

  int b_to_i(const std::array<uint8_t, S>& b) const;

 private:
  void generate_combinations();

  std::vector<std::array<uint8_t, S>> m_states;
};

template <size_t S, size_t N>
Combinations<S, N>::Combinations() {
  generate_combinations();
}

template <size_t S, size_t N>
const std::vector<std::array<uint8_t, S>>& Combinations<S, N>::states() const {
  return m_states;
}

template <size_t S, size_t N>
size_t Combinations<S, N>::size() const {
  return m_states.size();
}

template <size_t S, size_t N>
std::array<uint8_t, S> Combinations<S, N>::i_to_b(size_t i) const {
  ED_ASSERT(i < m_states.size());
  return m_states[i];
}

template <size_t S, size_t N>
int Combinations<S, N>::b_to_i(const std::array<uint8_t, S>& b) const {
  auto it = std::find(m_states.begin(), m_states.end(), b);
  ED_ASSERT(it != m_states.end());
  return std::distance(m_states.begin(), it);
}

template <size_t S, size_t N>
void Combinations<S, N>::generate_combinations() {
  std::array<uint8_t, S> v{};
  std::fill(v.begin() + S - N, v.end(), 1);
  do {
    m_states.push_back(v);
  } while (std::next_permutation(v.begin(), v.end()));
}
