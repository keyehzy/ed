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

// 1D Hubbard model in momentum space with S sites and F particles.
// t is nearest-neighbors hopping, u is local interaction and mu is chemical
// potential. kinect energy is diagonal with e(k)=-2*t*cos(k), while the
// interaction is non-diagonal. k is restrict to k=2*pi*n/S, with n = 0, +-1,
// +-2, ..., +- F/2.
template <size_t S, size_t F>
class MomentumSpace {
  static_assert(S >= F);

 public:
  MomentumSpace(double t, double u, double mu);

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
  void construct_sparse_representation(double k);

  double m_t;
  double m_u;
  double m_mu;
  Combinations<S, F> m_combs{};
  SparseMatrix<double> m_sparse_representation{};
};

template <size_t S, size_t F>
MomentumSpace<S, F>::MomentumSpace(double t, double u, double mu)
    : m_t(t), m_u(u), m_mu(mu) {
  ED_ASSERT(m_combs.size() == binomial(S, F));
  construct_sparse_representation(0);
}

template <size_t S, size_t F>
size_t MomentumSpace<S, F>::size() const {
  return m_combs.size();
}

template <size_t S, size_t F>
size_t MomentumSpace<S, F>::dim() const {
  return size() * size();
}

template <size_t S, size_t F>
double MomentumSpace<S, F>::t() const {
  return m_t;
}

template <size_t S, size_t F>
double MomentumSpace<S, F>::u() const {
  return m_u;
}

template <size_t S, size_t F>
double MomentumSpace<S, F>::mu() const {
  return m_mu;
}

template <size_t S, size_t F>
Eigen::MatrixXd MomentumSpace<S, F>::dense() const {
  Eigen::MatrixXd hmat = Eigen::MatrixXd::Zero(dim(), dim());
  for (const auto& [index, val] : m_sparse_representation.elements()) {
    hmat(index.i, index.j) += val;
  }
  return hmat;
}

template <size_t S, size_t F>
Eigen::SparseMatrix<double> MomentumSpace<S, F>::sparse() const {
  Eigen::SparseMatrix<double> hmat(dim(), dim());
  for (const auto& [index, val] : m_sparse_representation.elements()) {
    hmat.coeffRef(index.i, index.j) += val;
  }
  return hmat;
}

template <size_t S, size_t F>
void MomentumSpace<S, F>::construct_sparse_representation(double k) {
  // In this case, since we have momentum space representation, the kinect
  // energy is diagonal and has values of -2*t*cos(k), where k=2*pi*n/S with
  // n=0,+-1,+-2,...,F/S. Notice that k is defined mod 2*pi and -pi and +pi
  // corresponds to same state.
  //
  // On the other hand, the interacting part becomes a momentum preversing
  // scattering process, which is not diagonal in this basis.
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
