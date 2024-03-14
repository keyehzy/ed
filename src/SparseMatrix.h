#pragma once

template <typename T>
class SparseMatrix {
 public:
  struct Index {
    size_t i;
    size_t j;

    bool operator==(const Index& other) const {
      return i == other.i && j == other.j;
    }
  };

  struct IndexHasher {
    static constexpr size_t p = 2654435761;
    size_t operator()(const Index& index) const {
      size_t ih = std::hash<size_t>{}(index.i);
      size_t jh = std::hash<size_t>{}(index.j);
      return (p + ih) * p + jh;
    }
  };

  SparseMatrix() {}

  const T& operator()(size_t i, size_t j) const { return m_lut[Index{i, j}]; }

  T& operator()(size_t i, size_t j) { return m_lut[Index{i, j}]; }

  const std::unordered_map<Index, T, IndexHasher>& elements() const {
    return m_lut;
  }

 private:
  std::unordered_map<Index, T, IndexHasher> m_lut{};
};
