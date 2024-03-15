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

template <size_t S>
std::vector<Nibble<S>> generate_nibble_transitions(const Nibble<S>& nibble) {
  std::vector<Nibble<S>> result;
  size_t n = nibble.size();
  for (size_t i = 0; i < n; i++) {
    if (nibble[i] == 1) {
      if (nibble[(i - 1 + n) % n] == 0) {
        Nibble<S> left = nibble;
        left[i] = 0;
        left[(i - 1 + n) % n] = 1;
        result.push_back(left);
      }
      if (nibble[(i + 1) % n] == 0) {
        Nibble<S> right = nibble;
        right[i] = 0;
        right[(i + 1) % n] = 1;
        result.push_back(right);
      }
    }
  }
  return result;
}

template <size_t S>
std::vector<State<S>> nearest_neighbors_hoppings(const State<S>& state) {
  std::vector<State<S>> result;
  for (const Nibble<S>& nibble_up :
       generate_nibble_transitions(state.up_nibble())) {
    result.emplace_back(nibble_up, state.down_nibble());
  }
  for (const Nibble<S>& nibble_down :
       generate_nibble_transitions(state.down_nibble())) {
    result.emplace_back(state.up_nibble(), nibble_down);
  }
  return result;
}

// Hubbard model in a one-dimensional chain of S sites and F particles.
// t is nearest-neighbors hopping, u is local interaction and mu is chemical
// potential.
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
  // The interaction term is diagonal in this basis, counting the double
  // occupation. Similarly the chemical potential.
  //
  // The kinect energy term is not diagonal and moves the states around to
  // neighboring sites.
  for (const Nibble<S>& up_nibble : combinations().states()) {
    for (const Nibble<S>& down_nibble : combinations().states()) {
      State<S> s(up_nibble, down_nibble);
      size_t i = combinations().b_to_i(up_nibble) * size() +
                 combinations().b_to_i(down_nibble);
      m_sparse_representation(i, i) += m_u * s.double_occ() - m_mu * s.size();
      for (const State<S>& state : nearest_neighbors_hoppings(s)) {
        size_t j = combinations().b_to_i(state.up_nibble()) * size() +
                   combinations().b_to_i(state.down_nibble());
        m_sparse_representation(i, j) += -m_t;
      }
    }
  }
}
