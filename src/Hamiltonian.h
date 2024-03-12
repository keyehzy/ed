#pragma once

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>

#include "Assert.h"
#include "Binomial.h"
#include "Combinations.h"
#include "State.h"

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE

template <size_t S, size_t F>
class Hamiltonian {
 public:
  Hamiltonian(double t, double u, double mu) : m_t(t), m_u(u), m_mu(mu) {
    ED_ASSERT(m_combs.size() == binomial(S, F));
  }

  size_t size() const { return m_combs.size(); }

  const Combinations<S, F>& combinations() const { return m_combs; }

  double t() const { return m_t; }

  double u() const { return m_u; }

  double mu() const { return m_mu; }

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix() const {
    size_t dim = size() * size();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hmat(dim, dim);
    for (const Nibble<S>& up_nibble : combinations().states()) {
      for (const Nibble<S>& down_nibble : combinations().states()) {
        State<S> s(up_nibble, down_nibble);
        size_t i = combinations().b_to_i(up_nibble) * size() +
                   combinations().b_to_i(down_nibble);
        hmat(i, i) = m_u * s.double_occ() - m_mu * s.size();
        for (const State<S>& state : s.hopping()) {
          size_t j = combinations().b_to_i(state.up_nibble()) * size() +
                     combinations().b_to_i(state.down_nibble());
          hmat(i, j) = -m_t;
        }
      }
    }
    return hmat;
  }

 private:
  Combinations<S, F> m_combs{};
  double m_t;
  double m_u;
  double m_mu;
};
