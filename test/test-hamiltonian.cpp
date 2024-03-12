#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "Hamiltonian.h"

using testing::ElementsAre;
using testing::IsEmpty;
using testing::UnorderedElementsAre;

TEST(test_hamiltonian, simple) {
  Hamiltonian<2, 1> system(/*t=*/1, /*u=*/1, /*mu=*/0);
  auto matrix = system.matrix();
}
