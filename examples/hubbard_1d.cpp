#include "Hamiltonian.h"

std::vector<double> linspace(double min, double max, size_t n = 50) {
  double delta = (max - min) / (n - 1);
  std::vector<double> result;
  for (size_t i = 0; i < n; i++) {
    result.push_back(min + delta * i);
  }
  return result;
}

int main() {
  Hamiltonian<2, 1> system(/*t=*/1, /*u=*/1000, /*mu=*/0);
  auto matrix = system.dense();
  std::cout << matrix << std::endl;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
  es.compute(matrix);
  std::cout << es.eigenvalues().transpose() << std::endl;
  std::cout << '\n';

  auto gs = es.eigenvectors().col(0).cwiseAbs();

  for (int i = 0; i < gs.size(); i++) {
    std::cout << system.combinations().i_to_b(i / system.size()) << " "
              << system.combinations().i_to_b(i % system.size()) << " " << gs[i]
              << std::endl;
  }
  return 0;
}

int main1() {
  for (double uval : linspace(0, 100, 10)) {
    Hamiltonian<2, 1> system(/*t=*/1, /*u=*/uval, /*mu=*/0);
    auto matrix = system.dense();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
    es.compute(matrix);
    auto gs = es.eigenvectors().col(0).cwiseAbs().transpose();
    std::cout << uval << " " << gs << std::endl;
  }

  return 0;
}
