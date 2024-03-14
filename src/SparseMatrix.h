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
    size_t operator()(const Index& index) const {
      size_t prime = 2654435761;
      return (prime + std::hash<size_t>{}(index.i)) * prime +
             std::hash<size_t>{}(index.j);
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

