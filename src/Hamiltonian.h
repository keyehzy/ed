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
#include "SparseMatrix.h"
#include "State.h"

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE

template <size_t S, size_t F>
class Hubbard1D {
  static_assert(S >= F);

 public:
  Hubbard1D(double t, double u, double mu);

  size_t size() const;

  size_t dim() const;

  double t() const;

  double u() const;

  double mu() const;

  const Combinations<S, F>& combinations() const;

  const SparseMatrix<double>& sparse_representation() const;

  Eigen::MatrixXd dense() const;

  Eigen::SparseMatrix<double> sparse() const;

 private:
  void construct_sparse_representation();

  double m_t;
  double m_u;
  double m_mu;
  Combinations<S, F> m_combs{};
  SparseMatrix<double> m_sparse_representation{};
};

template <size_t S, size_t F>
Hubbard1D<S, F>::Hubbard1D(double t, double u, double mu)
    : m_t(t), m_u(u), m_mu(mu) {
  ED_ASSERT(m_combs.size() == binomial(S, F));
  construct_sparse_representation();
}

template <size_t S, size_t F>
size_t Hubbard1D<S, F>::size() const {
  return m_combs.size();
}

template <size_t S, size_t F>
size_t Hubbard1D<S, F>::dim() const {
  return size() * size();
}

template <size_t S, size_t F>
double Hubbard1D<S, F>::t() const {
  return m_t;
}

template <size_t S, size_t F>
double Hubbard1D<S, F>::u() const {
  return m_u;
}

template <size_t S, size_t F>
double Hubbard1D<S, F>::mu() const {
  return m_mu;
}

template <size_t S, size_t F>
const Combinations<S, F>& Hubbard1D<S, F>::combinations() const {
  return m_combs;
}

template <size_t S, size_t F>
const SparseMatrix<double>& Hubbard1D<S, F>::sparse_representation() const {
  return m_sparse_representation;
}

template <size_t S, size_t F>
Eigen::MatrixXd Hubbard1D<S, F>::dense() const {
  Eigen::MatrixXd hmat = Eigen::MatrixXd::Zero(dim(), dim());
  for (const auto& [index, val] : m_sparse_representation.elements()) {
    hmat(index.i, index.j) += val;
  }
  return hmat;
}

template <size_t S, size_t F>
Eigen::SparseMatrix<double> Hubbard1D<S, F>::sparse() const {
  Eigen::SparseMatrix<double> hmat(dim(), dim());
  for (const auto& [index, val] : m_sparse_representation.elements()) {
    hmat.coeffRef(index.i, index.j) += val;
  }
  return hmat;
}

template <size_t S, size_t F>
void Hubbard1D<S, F>::construct_sparse_representation() {
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
