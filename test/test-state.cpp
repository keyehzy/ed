#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "State.h"

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
