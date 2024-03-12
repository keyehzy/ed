#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "Hamiltonian.h"

using testing::ElementsAre;
using testing::IsEmpty;
using testing::UnorderedElementsAre;

TEST(test_hamiltonian, simple) { Hamiltonian<4, 2> system(1, 1, 0); }
