#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "Hamiltonian.h"

using testing::DoubleNear;
using testing::Pointwise;
using testing::UnorderedElementsAre;

TEST(test_hamiltonian, hubbard1d_groundstate) {
  {
    Hubbard1D<2, 1> system(/*t=*/1, /*u=*/0, /*mu=*/0);
    auto matrix = system.dense();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(matrix);
    auto gs = es.eigenvectors().col(0).cwiseAbs().transpose();
    EXPECT_THAT(gs, Pointwise(DoubleNear(1e-9), {0.5, 0.5, 0.5, 0.5}));
  }

  {
    Hubbard1D<2, 1> system(/*t=*/1, /*u=*/1e12, /*mu=*/0);
    auto matrix = system.dense();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(matrix);
    auto gs = es.eigenvectors().col(0).cwiseAbs().transpose();
    EXPECT_THAT(gs, Pointwise(DoubleNear(1e-9),
                              {0.0, 0.70710678118654, 0.70710678118654, 0.0}));
  }
}

TEST(test_hamiltonian, dense_vs_sparse) {
  Hubbard1D<2, 1> system(/*t=*/1, /*u=*/1, /*mu=*/0);

  // Dense
  auto dense = system.dense();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  es.compute(dense);

  // Sparse
  auto sparse = system.sparse();
  Spectra::SparseSymMatProd<double> op(sparse);
  Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, 1, 2);
  eigs.init();
  eigs.compute(Spectra::SortRule::SmallestAlge);

  double w_dense = es.eigenvalues()[0];
  double w_sparse = eigs.eigenvalues()[0];
  EXPECT_NEAR(w_dense, w_sparse, 1e-9);

  Eigen::VectorXd v_dense = es.eigenvectors().col(0);
  Eigen::VectorXd v_sparse = eigs.eigenvectors().col(0);
  EXPECT_THAT(v_dense, Pointwise(DoubleNear(1e-9), v_sparse));
}
