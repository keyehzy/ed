#pragma once

#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/SparseCore>
#include <unordered_map>

#include "Assert.h"
#include "Binomial.h"
#include "Combinations.h"
#include "State.h"

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE

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

  void insert(size_t i, size_t j, T val) { m_lut.emplace(i, j); }

  const T& operator()(size_t i, size_t j) const { return m_lut[Index{i, j}]; }

  T& operator()(size_t i, size_t j) { return m_lut[Index{i, j}]; }

  const std::unordered_map<Index, T, IndexHasher>& elements() const {
    return m_lut;
  }

 private:
  std::unordered_map<Index, T, IndexHasher> m_lut{};
};

template <size_t S, size_t F>
class Hubbard1D {
  static_assert(S >= F);

 public:
  Hubbard1D(double t, double u, double mu) : m_t(t), m_u(u), m_mu(mu) {
    ED_ASSERT(m_combs.size() == binomial(S, F));
    construct_sparse_representation();
  }

  size_t size() const { return m_combs.size(); }

  size_t dim() const { return size() * size(); }

  double t() const { return m_t; }

  double u() const { return m_u; }

  double mu() const { return m_mu; }

  const Combinations<S, F>& combinations() const { return m_combs; }

  const SparseMatrix<double>& sparse_representation() const {
    return m_sparse_representation;
  }

  Eigen::MatrixXd dense() const {
    Eigen::MatrixXd hmat = Eigen::MatrixXd::Zero(dim(), dim());
    for (const auto& [index, val] : m_sparse_representation.elements()) {
      hmat(index.i, index.j) += val;
    }
    return hmat;
  }

  Eigen::SparseMatrix<double> sparse() const {
    Eigen::SparseMatrix<double> hmat(dim(), dim());
    for (const auto& [index, val] : m_sparse_representation.elements()) {
      hmat.coeffRef(index.i, index.j) += val;
    }
    return hmat;
  }

 private:
  void construct_sparse_representation() {
    for (const Nibble<S>& up_nibble : combinations().states()) {
      for (const Nibble<S>& down_nibble : combinations().states()) {
        State<S> s(up_nibble, down_nibble);
        size_t i = combinations().b_to_i(up_nibble) * size() +
                   combinations().b_to_i(down_nibble);
        m_sparse_representation(i, i) += m_u * s.double_occ() - m_mu * s.size();
        for (const State<S>& state : s.hopping()) {
          size_t j = combinations().b_to_i(state.up_nibble()) * size() +
                     combinations().b_to_i(state.down_nibble());
          m_sparse_representation(i, j) += -m_t;
        }
      }
    }
  }

  double m_t;
  double m_u;
  double m_mu;
  Combinations<S, F> m_combs{};
  SparseMatrix<double> m_sparse_representation{};
};
