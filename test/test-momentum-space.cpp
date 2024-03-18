#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MomentumSpace.h"

using testing::DoubleNear;
using testing::IsEmpty;
using testing::Pointwise;
using testing::UnorderedElementsAre;

TEST(test_momentum_space, transitions) {
  MomentumSpaceInteractionGenerator<2> generator{};

  {
    Nibble<2> up{0, 0}, down{0, 0};
    State<2> state(up, down);
    size_t k = 0;
    EXPECT_THAT(generator.generate_transitions(k, state), IsEmpty());
  }

  {
    Nibble<2> up{1, 1}, down{0, 0};
    State<2> state(up, down);
    size_t k = 0;
    std::vector<State<2>> transitions =
        generator.generate_transitions(k, state);
    EXPECT_THAT(transitions, IsEmpty());
  }

  {
    Nibble<2> up{0, 0}, down{1, 1};
    State<2> state(up, down);
    size_t k = 0;
    std::vector<State<2>> transitions =
        generator.generate_transitions(k, state);
    EXPECT_THAT(transitions, IsEmpty());
  }

  {
    Nibble<2> up{1, 1}, down{1, 1};
    State<2> state(up, down);
    size_t k = 0;
    std::vector<State<2>> transitions =
        generator.generate_transitions(k, state);
    EXPECT_EQ(transitions.size(), 2);
    EXPECT_THAT(transitions, UnorderedElementsAre(
                                 State<2>(Nibble<2>{1, 1}, Nibble<2>{1, 1}),
                                 State<2>(Nibble<2>{1, 1}, Nibble<2>{1, 1})));
  }

  {
    Nibble<2> up{1, 0}, down{1, 0};
    State<2> state(up, down);
    size_t k = 0;
    std::vector<State<2>> transitions =
        generator.generate_transitions(k, state);
    EXPECT_EQ(transitions.size(), 2);
    EXPECT_THAT(transitions, UnorderedElementsAre(
                                 State<2>(Nibble<2>{1, 0}, Nibble<2>{1, 0}),
                                 State<2>(Nibble<2>{0, 1}, Nibble<2>{0, 1})));
  }

  {
    Nibble<2> up{0, 1}, down{1, 0};
    State<2> state(up, down);
    size_t k = 0;
    std::vector<State<2>> transitions =
        generator.generate_transitions(k, state);
    EXPECT_EQ(transitions.size(), 2);
    EXPECT_THAT(transitions, UnorderedElementsAre(
                                 State<2>(Nibble<2>{0, 1}, Nibble<2>{1, 0}),
                                 State<2>(Nibble<2>{1, 0}, Nibble<2>{0, 1})));
  }

  {
    Nibble<2> up{1, 0}, down{0, 1};
    State<2> state(up, down);
    size_t k = 0;
    std::vector<State<2>> transitions =
        generator.generate_transitions(k, state);
    EXPECT_THAT(transitions, IsEmpty());
  }
}
