#include <iostream>
#include <complex>
// #include <string>
#include <gtest/gtest.h>
#include "fact.hpp"
#include "erfc.hpp"
#include "prim.hpp"
#include "lcomb.hpp"
#include "hatom.hpp"

using namespace erfc_mori;
using namespace l2func;
using namespace fact;

TEST(first, first) {

  EXPECT_EQ(2, 1+1);
  
}
TEST(Factorial, Factorial) {

  EXPECT_ANY_THROW(Factorial(-1));
  EXPECT_EQ(1, Factorial(0));
  EXPECT_EQ(1, Factorial(1));
  EXPECT_EQ(2, Factorial(2));
  EXPECT_EQ(6, Factorial(3));
  EXPECT_EQ(24, Factorial(4));
  EXPECT_EQ(120, Factorial(5));
  EXPECT_EQ(720, Factorial(6));

  EXPECT_ANY_THROW(DoubleFactorial(-1));
  EXPECT_EQ(1,   DoubleFactorial(0));
  EXPECT_EQ(1,   DoubleFactorial(1));
  EXPECT_EQ(2,   DoubleFactorial(2));
  EXPECT_EQ(3,   DoubleFactorial(3));
  EXPECT_EQ(8,   DoubleFactorial(4));
  EXPECT_EQ(15,  DoubleFactorial(5));
  EXPECT_EQ(48,  DoubleFactorial(6));
  EXPECT_EQ(105, DoubleFactorial(7));
  
  
}
TEST(Erfc, real) {


  // this function is forbidden
  //  erfc_add_Eh_q<int>(1, 2);
  
  double y, x, expect;
  ErfcCalcData calc_data;

  x = 1.0;
  expect =0.157299207050285130658779364917390740703933002034;
  //erfc_d(x, y, calc_data);
  Erfc<double>(x, y, calc_data);
  EXPECT_DOUBLE_EQ(expect, y);
  EXPECT_TRUE(calc_data.convergence);
  
  x = 1.0;
  x = 1.0 / 100.0;
  expect = 0.98871658444415038308409047645193078905089904517;
  //  erfc_d(x, y, calc_data);
  Erfc<double>(x, y, calc_data);
  EXPECT_DOUBLE_EQ(expect, y);
  EXPECT_TRUE(calc_data.convergence);

  x = 3.0;
  expect = 0.0000220904969985854413727761295823203798477070873992;
  //  erfc_d(x, y, calc_data);
  // erfc(x, y, calc_data);
  Erfc<double>(x, y, calc_data);
  EXPECT_DOUBLE_EQ(expect, y);
  EXPECT_TRUE(calc_data.convergence);
  
}
TEST(Erfc, complex) {
  
  CD x, y, y_expect;
  ErfcCalcData calc_data;
  double eps = 10.0 * machine_eps();

  x.real() = 1;
  x.imag() = -1;

  y_expect.real() = -0.31615128169794764488027108024367;
  y_expect.imag() = +0.190453469237834686284108861969162;

  Erfc(x, y, calc_data);

  EXPECT_DOUBLE_EQ(y.real(), y_expect.real());
  EXPECT_NEAR(y.real(), y_expect.real(), eps);
  EXPECT_DOUBLE_EQ(y.imag(), y_expect.imag());
  EXPECT_NEAR(y.imag(), y_expect.imag(), eps);


  x.real() = 0.0157073173118206757532953533099;
  x.imag() = 0.9998766324816605986389071277312;
  
  y_expect.real() = 0.95184545524179913420177473658805;
  y_expect.imag() =-1.64929108965086517748934245403332;
  Erfc(x, y, calc_data);
  EXPECT_TRUE(calc_data.convergence);
  EXPECT_NEAR(     y.real(), y_expect.real(), eps);
  EXPECT_NEAR(     y.imag(), y_expect.imag(), eps);

  x.real() = 0.001564344650402308690101053194671668923139;
  x.imag() = 0.009876883405951377261900402476934372607584;
  y_expect.real() = 0.9982346553205423153337357292658472915601;
  y_expect.imag() =-0.0111452046101524188315708507537751407281;
  Erfc(x, y, calc_data);
  EXPECT_TRUE(calc_data.convergence);
  EXPECT_NEAR(y.real(), y_expect.real(), eps);
  EXPECT_NEAR(y.imag(), y_expect.imag(), eps);
}
TEST(Prim, Construct) {

  RSTO s1(2, 1.0);
  RGTO g1(2, 1.0);
  RSTO s2(3, 1.1, Normalized);

  EXPECT_EQ(2, s1.n());
  EXPECT_EQ(2, g1.n());
  //  EXPECT_EQ(3, s2.n());

  CSTO s3(3, CD(1.0, -0.2));
  CGTO g3(3, CD(1.0, -0.2));
  EXPECT_EQ(3, s3.n());
  EXPECT_EQ(3, g3.n());
}
TEST(Prim, OrbitalExp) {
  RSTO s1(2, 1.1);
  s1.set_z(1.3);
  EXPECT_DOUBLE_EQ(1.3, s1.z());
}
TEST(Prim, CIP) {
  RSTO s1(2.5, 2, 1.1);
  CSTO s2(1.2, 3, CD(0.4, 0.2));
  RGTO g1(0.3, 1, 1.2);
  CGTO g2(0.4, 4, CD(0.1, -0.1));
  double eps = 10.0 * machine_eps();

  //  EXPECT_NEAR(2.9105687018397886, CIP(s1, s1), eps);
  //  EXPECT_NEAR(0.0764996359892135, CIP(s1, g1), eps);
  
  EXPECT_DOUBLE_EQ(2.9105687018397886, CIP(s1, s1));
  
  EXPECT_DOUBLE_EQ(20.98619989895233, CIP(CSTO(s1), s2).real());
  EXPECT_DOUBLE_EQ(20.98619989895233, CIP(s2, CSTO(s1)).real());
  EXPECT_DOUBLE_EQ(-21.40636768181864, CIP(CSTO(s1), s2).imag());
  EXPECT_DOUBLE_EQ(-21.40636768181864, CIP(s2, CSTO(s1)).imag());
  
  EXPECT_NEAR(0.0764996359892135, CIP(s1, g1), eps);
  EXPECT_NEAR(0.0764996359892135, CIP(g1, s1), eps);

  CSTO c_s1(s1);
  eps *= 100000;
  EXPECT_NEAR(5.562595882704702, CIP(c_s1, g2).real(),      eps);
  EXPECT_NEAR(5.562595882704702, CIP(g2, CSTO(s1)).real(),  eps);
  EXPECT_NEAR(+22.587241177071004, CIP(CSTO(s1), g2).imag(),eps);
  EXPECT_NEAR(+22.587241177071004, CIP(g2, CSTO(s1)).imag(),eps);

  EXPECT_NEAR(0.05270913901892936, CIP(CGTO(g1), g2).real(), eps);
  EXPECT_NEAR(0.05270913901892936, CIP(g2, CGTO(g1)).real(), eps);
  EXPECT_NEAR(0.012359047425198447, CIP(CGTO(g1), g2).imag(), eps);
  EXPECT_NEAR(0.012359047425198447, CIP(g2, CGTO(g1)).imag(), eps);
  
}
TEST(Prim, Normalized) {
  
  RSTO n_s1(2, 1.1, Normalized);
  RSTO s1(1.0, 2, 1.1);
  RSTO n_s2(3, 1.2, Normalized);
  RSTO s2(1.0, 3, 1.2);
  
  EXPECT_DOUBLE_EQ(1.0, CIP(n_s1, n_s1));

  EXPECT_NEAR( CIP(n_s1, Op(OpDDr<RSTO>(), n_s2)),
	       CIP(s1, Op(OpDDr<RSTO>(), s2)) /
	       sqrt(CIP(s1, s1) * CIP(s2, s2)),
	       0.000000000001);

  RGTO n_g1(2, 1.1, Normalized);
  RGTO g1(1.0, 2, 1.1);
  RGTO n_g2(3, 1.2, Normalized);
  RGTO g2(1.0, 3, 1.2);
  
  EXPECT_DOUBLE_EQ(1.0, CIP(n_g1, n_g1));

  EXPECT_NEAR( CIP(n_g1, Op(OpDDr<RGTO>(), n_g2)),
	       CIP(g1,   Op(OpDDr<RGTO>(), g2)) /
	       sqrt(CIP(g1, g1) * CIP(g2, g2)),
	       0.000000000001);  
}
TEST(Prim, power) {

  EXPECT_EQ(1, CSTO::exp_power);
  
}
TEST(Prim, DBasis) {
  
  RSTO d_s1 = DBasis<1, RSTO>(1.0, 3, 0.2);
  EXPECT_DOUBLE_EQ(-1.0, d_s1.c());
  EXPECT_EQ(4, d_s1.n());
  EXPECT_EQ(0.2, d_s1.z());
  
}
TEST(Prim, AtX) {

  RSTO s1(1.1, 2, 0.2);
  double x0(3.0);
  EXPECT_DOUBLE_EQ(1.1 * x0 * x0 * exp(-0.2 * x0), AtX(x0, s1));
  
}
TEST(LinearComb, Construct) {

  LinearComb<RSTO> rstos;
  //  rstos.add(1.3, RSTO(2, 1.5));
  //  rstos.add(1.3 * RSTO(2, 1.5));
  rstos += 1.3 * RSTO(2, 1.5);
  EXPECT_DOUBLE_EQ(1.3, rstos.coef_i(0));

  rstos[0] = 1.1 * RSTO(3, 1.1);
  EXPECT_DOUBLE_EQ(1.1, rstos.coef_i(0));

  LinearComb<RSTO> f;
  f += rstos; 
  f += RSTO(1.0, 1, 1.1);
  EXPECT_EQ(2, f.size());
  
}
TEST(LinearComb, CIP) {

  CSTO s1(1.1, 2, 3.1);
  CSTO s2(1.2, 3, 0.1);
  CSTO s3(1.2, 3, 0.1);

  LinearComb<CSTO> f1;
  LinearComb<CSTO> f2;
  f1 += 1.2 * s1;
  f1 += 1.3 * s2;
  f2 += 0.3 * s3;

  CD expe = 0.3 * 1.2 * CIP(s1, s3) + 0.3 * 1.3 * CIP(s2, s3);
  CD calc = CIP(f1, f2);
  EXPECT_DOUBLE_EQ(expe.real(), calc.real());
  EXPECT_DOUBLE_EQ(expe.imag(), calc.imag());

  EXPECT_DOUBLE_EQ(CIP(s1, f2).real(), 0.3 * CIP(s1, s3).real());
  EXPECT_DOUBLE_EQ(CIP(f2, s1).real(), 0.3 * CIP(s1, s3).real());

}
TEST(LinearComb, op) {

  LinearComb<RSTO> f1;
  f1 += 1.2 * RSTO(1.1, 2, 0.2);
  f1 += 1.3 * RSTO(0.9, 3, 0.3);

  LinearComb<RSTO> f2 = Op(OpRm<RSTO>(2), f1);
  EXPECT_EQ(2, f2.size());
  EXPECT_EQ(4, f2.prim_i(0).n());
  EXPECT_EQ(5, f2.prim_i(1).n());
  EXPECT_DOUBLE_EQ(1.3, f2.coef_i(1));
  
  LinearComb<RSTO> f3 = Op(OpRm<RSTO>(3), RSTO(1.0, 2, 1.0));
  EXPECT_EQ(1, f3.size());
  EXPECT_EQ(5, f3.prim_i(0).n());

  LinearComb<RSTO> df1 = Op(OpDDr<RSTO>(),
			    RSTO(1.1, 2, 1.2));
  EXPECT_EQ(2, df1.size());

  LinearComb<RSTO> df2 = Op(OpDDr2<RSTO>(),
			    RSTO(1.1, 2, 1.2));
  EXPECT_EQ(4, df2.size());  
  
}
TEST(LinearComb, AtX) {

  RSTO s1(1.1, 1, 0.2);
  RSTO s2(1.2, 1, 0.3);
  LinearComb<RSTO> sto;
  sto += 0.3 * s1;
  sto += 0.4 * s2;
  
  double x0(3.0);
  EXPECT_DOUBLE_EQ(x0 * (0.3 * 1.1 * exp(-0.2 * x0) +
			 0.4 * 1.2 * exp(-0.3*x0)),
		   AtX(x0, sto));
  
}
TEST(LinearComb, dbasis) {

  RSTO s1(1.2, 2, 0.3);
  LinearComb<RSTO> d_s1  = D1Normalized(s1);
  LinearComb<RSTO> dd_s1 = D2Normalized(s1);

  double x0 = 3.3;
  double dx = 0.01;
  EXPECT_NEAR( (AtX(x0+dx, s1) + AtX(x0-dx, s1) - 2 * AtX(x0, s1)) / dx,
	       AtX(x0, d_s1),  0.00001);
  EXPECT_NEAR( (AtX(x0 + dx, s1) + AtX(x0 - dx, s1))/(dx * dx),
	       AtX(x0, dd_s1), 0.00001);  
}
TEST(HAtom, EigenState) {

  LinearComb<CSTO> f00 = HAtomEigenState<CD>(1, 0);
  LinearComb<CSTO> f10 = HAtomEigenState<CD>(2, 0);
  LinearComb<CSTO> f11 = HAtomEigenState<CD>(2, 1);
  EXPECT_EQ(1, f00.size());

  EXPECT_EQ(0.0, abs(CIP(f00, f10)));
  EXPECT_EQ(0.0, abs(CIP(f00, f10)));

  EXPECT_DOUBLE_EQ(1.0, CIP(f00, f00).real());
  EXPECT_DOUBLE_EQ(1.0, CIP(f10, f10).real());
  EXPECT_DOUBLE_EQ(1.0, CIP(f11, f11).real());

  double eps = machine_eps() * 10;
  
  EXPECT_NEAR(HAtomEigenEnergy(1),
	      CIP(f00, Op(OpHAtomRadial<CSTO>(1.0, 0), f00)).real(), eps);
  EXPECT_NEAR(0.0,
	      CIP(f00, Op(OpHAtomRadial<CSTO>(1.0, 0), f10)).real(), eps);
  EXPECT_NEAR(HAtomEigenEnergy(2),
	      CIP(f10, Op(OpHAtomRadial<CSTO>(1.0, 0), f10)).real(), eps);
  EXPECT_NEAR(HAtomEigenEnergy(2),
	      CIP(f11, Op(OpHAtomRadial<CSTO>(1.0, 1), f11)).real(), eps);
}
TEST(HAtom, DipoleInitLength) {

  LinearComb<CSTO> mu_phi_00 = HAtomDipoleInitL<CD>(1, 0, 1);
  EXPECT_EQ(1, mu_phi_00.size());
  
}



