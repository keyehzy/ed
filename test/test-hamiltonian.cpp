#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "Hamiltonian.h"

using testing::DoubleNear;
using testing::Pointwise;

TEST(test_hamiltonian, dense_vs_sparse) {
  Hamiltonian<2, 1> system(/*t=*/1, /*u=*/1, /*mu=*/0);

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
