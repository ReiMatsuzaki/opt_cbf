#ifndef ADD_GTEST_TEMPLATE_H
#define ADD_GTEST_TEMPLATE_H

#include <complex>
#include <gtest/gtest.h>

::testing::AssertionResult AssertComplexEq(const char* a_expr,
					   const char* b_expr,
					   std::complex<double> a,
					   std::complex<double> b) {

  typedef ::testing::internal::FloatingPoint<double> FP;
  const FP a_r(a.real()), a_i(a.imag()), b_r(b.real()), b_i(b.imag());

  if(a_r.AlmostEquals(b_r) && a_i.AlmostEquals(b_i)) 
    return ::testing::AssertionSuccess();

  return ::testing::AssertionFailure()
    << a_expr << " and " << b_expr << " are not near." << std::endl
    << a_expr << " : " << a << std::endl 
    << b_expr << " : " << b << std::endl
    << "|" << a_expr << "-" << b_expr<< "| : " << abs(a-b) << std::endl;

}

::testing::AssertionResult AssertComplexNear(const char* a_expr,
					     const char* b_expr,
					     const char* eps_expr,
					     std::complex<double> a,
					     std::complex<double> b,
					     double eps) {


  if(abs(a-b) < eps)
    return ::testing::AssertionSuccess();

  return ::testing::AssertionFailure()
    << a_expr << " and " << b_expr << " are not near." << std::endl
    << a_expr << " : " << a << std::endl 
    << b_expr << " : " << b << std::endl
    << eps_expr << " : " << eps << std::endl
    << "|" << a_expr << "-" << b_expr<< "| : " << abs(a-b) << std::endl;

}

#define EXPECT_C_EQ(a, b) EXPECT_PRED_FORMAT2(AssertComplexEq, a, b)
#define EXPECT_C_NEAR(a, b, c) EXPECT_PRED_FORMAT3(AssertComplexNear, (a), (b), (c))


#endif
