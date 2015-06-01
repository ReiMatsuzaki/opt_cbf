#include <gtest/gtest.h>
#include "keys_values.hpp"

TEST(first, first) {
  EXPECT_EQ(2, 1 + 1);
}
TEST(TestKeysValues, Construct) {
  KeysValues keys_values;

  keys_values.Add("a", 1.1);
  keys_values.Add("a", 1.2);
  EXPECT_TRUE(keys_values.ExistKey("a"));
  EXPECT_EQ(2, keys_values.Count("a"));
  EXPECT_DOUBLE_EQ(1.1, keys_values.Get<double>("a", 0));
  EXPECT_NO_THROW(keys_values.Get<double>("a", 1));
  EXPECT_DOUBLE_EQ(1.2, keys_values.Get<double>("a", 1));
  EXPECT_ANY_THROW(keys_values.Get<double>("a", 2));

  keys_values.Add("b", 1);
  EXPECT_ANY_THROW(keys_values.Add("b", 1.2));

  keys_values.Add("c", make_tuple(1, 1.2));
  keys_values.Add("c", make_tuple(2, 1.3));
  EXPECT_ANY_THROW(keys_values.Add("c", 3));
  
}
