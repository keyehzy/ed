#pragma once

template <size_t S>
using Nibble = std::array<uint8_t, S>;

template <size_t S>
std::ostream& operator<<(std::ostream& out, const Nibble<S>& nibble) {
  for (uint8_t x : nibble) {
    out << (x ? 1 : 0);
  }
  return out;
}

template <size_t S>
class State {
 public:
  State(const Nibble<S>& up_nibble, const Nibble<S>& down_nibble)
      : m_up_nibble(up_nibble), m_down_nibble(down_nibble) {}

  const Nibble<S>& up_nibble() const { return m_up_nibble; }

  const Nibble<S>& down_nibble() const { return m_down_nibble; }

  size_t count_spin_down() const {
    size_t cnt = 0;
    for (size_t i = 0; i < size(); i++) {
      cnt += m_down_nibble[i] ? 1 : 0;
    }
    return cnt;
  }

  size_t count_spin_up() const {
    size_t cnt = 0;
    for (size_t i = 0; i < size(); i++) {
      cnt += m_up_nibble[i] ? 1 : 0;
    }
    return cnt;
  }

  size_t double_occ() const {
    size_t cnt = 0;
    for (size_t i = 0; i < size(); i++) {
      cnt += m_down_nibble[i] && m_up_nibble[i] ? 1 : 0;
    }
    return cnt;
  }

  size_t size() const { return S; }

  bool operator==(const State& other) const {
    return m_down_nibble == other.down_nibble() &&
           m_up_nibble == other.up_nibble();
  }

  friend std::ostream& operator<<(std::ostream& out, const State& state) {
    return out << state.up_nibble() << state.down_nibble();
  }

 private:
  Nibble<S> m_up_nibble;
  Nibble<S> m_down_nibble;
};
