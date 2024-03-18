#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/SparseCore>
#include <numbers>
#include <unordered_map>

#include "Assert.h"
#include "Binomial.h"
#include "Combinations.h"
#include "SparseMatrix.h"
#include "State.h"

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE

template <size_t S>
struct MomentumSpaceInteractionGenerator {
  constexpr bool can_move_nibble(Nibble<S> nibble, size_t i, size_t j,
                                 size_t n) {
    return nibble[i] == 1 && (nibble[(j + n) % n] == 0 ||
                              (nibble[(j + n) % n] == 1 && i == (j + n) % n));
  }

  constexpr std::vector<State<S>> generate_transitions(size_t k,
                                                       const State<S>& state) {
    Nibble<S> up_nibble = state.up_nibble();
    Nibble<S> down_nibble = state.down_nibble();
    ED_ASSERT(up_nibble.size() == down_nibble.size());
    std::vector<State<S>> result;
    size_t n = up_nibble.size();
    for (size_t p = 0; p < n; p++) {
      for (size_t q = 0; q < n; q++) {
        Nibble<S> down = down_nibble;
        Nibble<S> up = up_nibble;
        bool to_push = false;
        if (can_move_nibble(up_nibble, p, p + q, n) &&
            can_move_nibble(down_nibble, k, k - q, n)) {
          up[p] = 0;
          up[(p + q + n) % n] = 1;
          down[k] = 0;
          down[(k - q + n) % n] = 1;
          result.emplace_back(up, down);
        }
      }
    }
    return result;
  }
};

// 1D Hubbard model in momentum space with S sites and F particles.
// t is nearest-neighbors hopping, u is local interaction and mu is chemical
// potential. kinect energy is diagonal with e(k)=-2*t*cos(k), while the
// interaction is non-diagonal. k is restrict to k=2*pi*n/S, with n = 0, +-1,
// +-2, ..., +- F/2.
template <size_t S, size_t F>
class MomentumSpace {
  static_assert(S >= F);

 public:
  MomentumSpace(double t, double u, double mu) : m_t(t), m_u(u), m_mu(mu) {
    ED_ASSERT(m_combs.size() == binomial(S, F));
    construct_sparse_representation(0);
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

  Eigen::MatrixXd dense() const;

  Eigen::SparseMatrix<double> sparse() const;

 private:
  void construct_sparse_representation(size_t n);

  double m_t;
  double m_u;
  double m_mu;
  Combinations<S, F> m_combs{};
  SparseMatrix<double> m_sparse_representation{};
  MomentumSpaceInteractionGenerator<S> m_interaction_generator{};
};

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

static inline double epsk(double t, size_t n, size_t N) {
  double k =
      2.0 * std::numbers::pi * static_cast<double>(n) / static_cast<double>(N);
  return -2.0 * t * std::cos(k);
}

static inline double local_u(double u, size_t N) {
  return u / static_cast<double>(N);
}

template <size_t S, size_t F>
void MomentumSpace<S, F>::construct_sparse_representation(size_t n) {
  // In this case, since we have momentum space representation, the kinect
  // energy is diagonal and has values of -2*t*cos(k), where k=2*pi*n/S with
  // n=0,+-1,+-2,...,F/2. There are in total 2*F/2+1=F+1 states. Notice that k
  // is defined mod 2*pi and -pi and +pi corresponds to same state.
  //
  // On the other hand, the interacting part becomes a momentum preversing
  // scattering process, which is not diagonal in this basis.
  for (const Nibble<S>& up_nibble : combinations().states()) {
    for (const Nibble<S>& down_nibble : combinations().states()) {
      State<S> s(up_nibble, down_nibble);
      size_t i = combinations().b_to_i(up_nibble) * size() +
                 combinations().b_to_i(down_nibble);
      m_sparse_representation(i, i) += epsk(m_t, n, size());
      for (const State<S>& state :
           m_interaction_generator.generate_transitions(n, s)) {
        size_t j = combinations().b_to_i(state.up_nibble()) * size() +
                   combinations().b_to_i(state.down_nibble());
        m_sparse_representation(i, j) += local_u(m_u, size());
      }
    }
  }
}
