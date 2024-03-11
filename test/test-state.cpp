#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "State.h"

using testing::ElementsAre;
using testing::IsEmpty;
using testing::UnorderedElementsAre;

TEST(test_state, simple) {
  Nibble<4> upnibble{1, 1, 0, 0};
  Nibble<4> downnibble{1, 0, 1, 0};
  State s(upnibble, downnibble);
  EXPECT_EQ(s.up_nibble(), upnibble);
  EXPECT_EQ(s.down_nibble(), downnibble);
  EXPECT_EQ(s.count_spin_up(), 2);
  EXPECT_EQ(s.count_spin_down(), 2);
  EXPECT_EQ(s.double_occ(), 1);
  EXPECT_EQ(s.size(), 4);
}

TEST(test_state, transitions) {
  {
    Nibble<4> upnibble{0, 0, 0, 0};
    Nibble<4> downnibble{0, 0, 0, 0};
    State s(upnibble, downnibble);
    EXPECT_THAT(s.hopping(), IsEmpty());
  }

  {
    Nibble<4> upnibble{1, 1, 1, 1};
    Nibble<4> downnibble{1, 1, 1, 1};
    State s(upnibble, downnibble);
    EXPECT_THAT(s.hopping(), IsEmpty());
  }

  {
    Nibble<4> up_nibble{1, 1, 0, 0};
    Nibble<4> down_nibble{1, 0, 1, 0};
    State s(up_nibble, down_nibble);
    EXPECT_THAT(s.hopping(), UnorderedElementsAre(
                                 State(up_nibble, Nibble<4>{0, 0, 1, 1}),
                                 State(up_nibble, Nibble<4>{0, 1, 1, 0}),
                                 State(up_nibble, Nibble<4>{1, 1, 0, 0}),
                                 State(up_nibble, Nibble<4>{1, 0, 0, 1}),
                                 State(Nibble<4>{0, 1, 0, 1}, down_nibble),
                                 State(Nibble<4>{1, 0, 1, 0}, down_nibble)));
  }
}
