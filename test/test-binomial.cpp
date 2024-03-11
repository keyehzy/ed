#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "Binomial.h"

using testing::ElementsAre;
using testing::IsEmpty;
using testing::UnorderedElementsAre;

// https://oeis.org/A000217
TEST(test_binomial, triangle_numbers) {
  EXPECT_EQ(binomial(2, 2), 1);
  EXPECT_EQ(binomial(3, 2), 3);
  EXPECT_EQ(binomial(4, 2), 6);
  EXPECT_EQ(binomial(5, 2), 10);
  EXPECT_EQ(binomial(6, 2), 15);
  EXPECT_EQ(binomial(7, 2), 21);
  EXPECT_EQ(binomial(8, 2), 28);
  EXPECT_EQ(binomial(9, 2), 36);
  EXPECT_EQ(binomial(10, 2), 45);
}
