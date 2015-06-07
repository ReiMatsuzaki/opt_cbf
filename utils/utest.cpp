#ifndef GTEST_HAS_TR1_TUPLE
#define GTEST_HAS_TR1_TUPLE 0

#include <complex>
#include <gtest/gtest.h>
#include "keys_values.hpp"
#include "timer.hpp"

TEST(TestConvertData, Atomic) {

  EXPECT_EQ(1, ConvertData<int>("1"));
  EXPECT_DOUBLE_EQ(1.1, ConvertData<double>("1.1"));
  EXPECT_DOUBLE_EQ(1.1,
		   ConvertData<CD>("(1.1,1.3)").real());
  
}
TEST(TestConvertData, Tuple) {

  EXPECT_EQ(1, 
	    get<1>(ConvertData<int,int>("34,1", ",")));
  EXPECT_DOUBLE_EQ(1.3,
		   get<1>(ConvertData<int,double>("1,1.3",",")));

  EXPECT_DOUBLE_EQ(1.3,
		   get<2>(ConvertData<int,double,double>
			  ("1,2,1.3", ",")));

  EXPECT_DOUBLE_EQ(1.3,
		   get<1>(ConvertData<int,CD,double>
			  ("2  (1.3,0.2) 5", " ")).real());

  EXPECT_DOUBLE_EQ(1.3,
		   get<1>(ConvertData<int,CD,double>
			  ("2  (1.3,0.2) 5", " ")).real());
}
TEST(TestKeysValues, Add) {
  
  KeysValues keys_values(":", " ");
  
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
TEST(TestKeysValues, AddAtom) {

  KeysValues keys_values(":", " ");

  keys_values.AddAtomConverting("a", "1.2");
  keys_values.AddAtomConverting("a", "2.2");
  keys_values.AddAtomConverting("a", "3.2");

  EXPECT_DOUBLE_EQ(2.2, keys_values.Get<double>("a", 1));
  EXPECT_ANY_THROW(keys_values.Get<int>("a", 2));
  EXPECT_ANY_THROW(keys_values.Get<double>("a", 3));

  keys_values.AddAtomConverting("bb", "3");
  EXPECT_ANY_THROW(keys_values.AddAtomConverting("bb", "1.1"));
}
TEST(TestKeysValues, ConvertData) {

  KeysValues keys_values(":", " ");

  keys_values.Add<string>("a", "1.2");
  keys_values.Add<string>("a", "2.2");
  keys_values.Add<string>("a", "3.2");

  keys_values.Add<string>("b", "(1.0,2.0)");
  keys_values.Add<string>("b", "(1.1,2.0)");
  keys_values.Add<string>("b", "2.0");

  keys_values.Add<string>("cc", "abc 1.03 (1.1,2.0)");
  keys_values.Add<string>("cc", "a 1.1 (1.0,8.0)");

  keys_values.ConvertValues<double>("a");
  keys_values.ConvertValues<CD>("b");
  keys_values.ConvertValues<string,double,CD>("cc");

  EXPECT_DOUBLE_EQ(2.2, keys_values.Get<double>("a", 1));
  EXPECT_NEAR(1.1, keys_values.Get<CD>("b", 1).real(),
	      0.0000001);
  EXPECT_NEAR(2.0, keys_values.Get<CD>("b", 2).real(),
	      0.0000001);
  EXPECT_DOUBLE_EQ(1.1,
		   get<1>(keys_values.Get<tuple<string,double,CD> >("cc", 1)));
}
TEST(TestKeysValues, SetIfNull) {

  KeysValues keys_values(":", " ");

  keys_values.Add<string>("a", "1.2");
  keys_values.Add<string>("a", "3.2");
  keys_values.SetIfNull<int>("b", 21);

  EXPECT_EQ(21, keys_values.Get<int>("b"));
}
TEST(TestKeysValues, ReadLine) {

  KeysValues keys_values(":", " ");

  keys_values.ReadLine("aa: 1.2  1.3  bbb  ");
  keys_values.ReadLine("  aa :  0.2    2.3   bcc  ");

  keys_values.ConvertValues<double, double, string>("aa");
  typedef tuple<double, double, string> TT;
  TT t = keys_values.Get<TT>("aa", 0);
  EXPECT_NEAR(1.2, get<0>(t),
	      0.0000001);
  EXPECT_EQ("bcc", get<2>(keys_values.Get<TT>("aa", 1)));
  
}
TEST(TestKeysValues, Read) {

  ifstream ifs("sample.in");
  KeysValues keys_values(":", " ");
  keys_values.Read(ifs);
  keys_values.ConvertValues<double,CD,string>("key1");
  keys_values.ConvertValues<double>("key2");
  keys_values.ConvertValues<string>("key3");

  typedef tuple<double, CD, string> TT;
  EXPECT_DOUBLE_EQ(2.1, get<1>(keys_values.Get<TT>("key1", 0)).imag());
  EXPECT_DOUBLE_EQ(8.3, keys_values.Get<double>("key2", 0));
  EXPECT_DOUBLE_EQ(1.0, keys_values.Get<double>("key2", 2));

  
}
TEST(TestTimer, Construct) {
  Timer timer;
  timer.Start("a");
  EXPECT_ANY_THROW(timer.End("b"));
  timer.End("a");
  timer.Display();
}
#endif
