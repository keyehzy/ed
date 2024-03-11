#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "Combinations.h"

using testing::ElementsAre;
using testing::IsEmpty;
using testing::UnorderedElementsAre;

TEST(test_combinations, corner_cases) {
  {
    Combinations<0, 0> combs;
    EXPECT_EQ(combs.size(), 1);
    EXPECT_THAT(combs.states(), ElementsAre(std::array<uint8_t, 0>{}));
  }

  {
    Combinations<1, 0> combs;
    EXPECT_EQ(combs.size(), 1);
    EXPECT_THAT(combs.states(), ElementsAre(std::array<uint8_t, 1>{0}));
  }

  {
    Combinations<1, 1> combs;
    EXPECT_EQ(combs.size(), 1);
    EXPECT_THAT(combs.states(), ElementsAre(std::array<uint8_t, 1>{1}));
  }

  {
    Combinations<2, 1> combs;
    EXPECT_EQ(combs.size(), 2);
    EXPECT_THAT(combs.states(),
                UnorderedElementsAre(std::array<uint8_t, 2>{1, 0},
                                     std::array<uint8_t, 2>{0, 1}));
  }

  {
    Combinations<2, 2> combs;
    EXPECT_EQ(combs.size(), 1);
    EXPECT_THAT(combs.states(),
                UnorderedElementsAre(std::array<uint8_t, 2>{1, 1}));
  }
}

TEST(test_combinations, elements) {
  {
    Combinations<4, 2> combs;
    EXPECT_EQ(combs.size(), 6);
    EXPECT_THAT(combs.states(),
                UnorderedElementsAre(std::array<uint8_t, 4>{0, 0, 1, 1},
                                     std::array<uint8_t, 4>{0, 1, 0, 1},
                                     std::array<uint8_t, 4>{1, 0, 0, 1},
                                     std::array<uint8_t, 4>{0, 1, 1, 0},
                                     std::array<uint8_t, 4>{1, 0, 1, 0},
                                     std::array<uint8_t, 4>{1, 1, 0, 0}));
  }
}

TEST(test_combinations, indices) {
  {
    Combinations<4, 2> combs;
    for (const std::array<uint8_t, 4>& state : combs.states()) {
      EXPECT_EQ(state, combs.i_to_b(combs.b_to_i(state)));
    }
  }
}
