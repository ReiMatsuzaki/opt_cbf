#include <iostream>
#include <complex>
#include <string>
#include <vector>
#include <boost/function.hpp>
#include <gtest/gtest.h>
#include <Eigen/Core>
#include <l2func.hpp>
#include <keys_values.hpp>
#include "l_algebra.hpp"
#include "restrict.hpp"
#include "opt.hpp"
#include "driv.hpp"
#include "opt_cbf.hpp"
#include "factory.hpp"

using namespace Eigen;
using namespace opt_cbf_h;
using namespace l2func;
using namespace std;

TEST(TEST_Create, CheckExceptions) {

  KeysValues kv(":", " ");

  // check throw exception when there does not exist
  // any basis functions data.
  EXPECT_ANY_THROW(CreateFactory(kv));

  kv.Add("opt_basis", make_tuple(1, CD(1.2, 0.0)));
  kv.Add("opt_basis", make_tuple(2, CD(3.4, 0.0)));
  kv.Add("opt_et_basis", make_tuple(1,5,
				    CD(1.2,0.0),CD(1.3,0.0)));

  // check throw exception when both opt_basis and opt_et_basis
  // exists simultouneously.
  EXPECT_ANY_THROW(CreateFactory(kv));
}
class TEST_MonoBasis : public ::testing::Test {
public:
  KeysValues kv;
  TEST_MonoBasis() : kv(":", " ") {}
  virtual void SetUp() {
    kv.Add<string>("channel", "1s->kp");
    kv.Add<string>("dipole", "length");
    kv.Add<double>("energy", 0.5);
    kv.Add<string>("basis_type", "GTO");
    kv.Add("opt_basis", make_tuple(1, CD(1.2, 0.0)));
    kv.Add("opt_basis", make_tuple(2, CD(3.4, 0.0)));
    kv.Add("opt_basis", make_tuple(2, CD(1.4, 0.0)));
    kv.Add("max_iter", 100);
    kv.Add("eps", 0.001);
  }

};
TEST_F(TEST_MonoBasis, Construct) {

  EXPECT_NO_THROW(CreateFactory(kv));
  EXPECT_TRUE(CreateFactory(kv) != NULL);

}
TEST_F(TEST_MonoBasis, BasisSet) {
  
  IFactory* factory = CreateFactory(kv);
  EXPECT_TRUE(factory != NULL);
  vector<CGTO>* gtos(NULL);
  EXPECT_ANY_THROW(factory->STOSet());

  gtos = factory->GTOSet();

  EXPECT_TRUE(gtos != NULL);
  
  EXPECT_EQ(3, gtos->size());
  EXPECT_EQ(2, (*gtos)[1].n());
  EXPECT_DOUBLE_EQ(1.4, (*gtos)[2].z().real());    
    
}
TEST_F(TEST_MonoBasis, OptTarget) {

  IFactory* factory = CreateFactory(kv);
  EXPECT_TRUE(factory != NULL);

  HAtomPI<l2func::CGTO>* hatom = NULL;
  
  hatom = factory->HAtomPiGTO();

  EXPECT_TRUE(hatom != NULL);

}
TEST_F(TEST_MonoBasis, Optimizer) {
  
  IFactory* factory = CreateFactory(kv);
  IOptimizer<CD>* opt(NULL);
  opt = factory->Optimizer();

  EXPECT_TRUE(opt != NULL);
  EXPECT_TRUE(typeid(*opt) == typeid(OptimizerNewton<CD>));
  EXPECT_EQ(100, opt->max_iter());
  EXPECT_DOUBLE_EQ(0.001, opt->eps());

}
class TEST_EvenTemp: public ::testing::Test {
public:
  // -------- Field ----------
  KeysValues kv;

  // -------- Method ---------
  TEST_EvenTemp() : kv(":", " ") {}
  virtual void SetUp() {
    kv.Add<string>("channel", "1s->kp");
    kv.Add<string>("dipole", "length");
    kv.Add<double>("energy", 0.5);
    kv.Add<string>("basis_type", "STO");
    kv.Add("opt_et_basis", make_tuple(1, 5, CD(1.2, 0.1), CD(1.3, 0.0)));
    kv.Add("opt_et_basis", make_tuple(2, 4, CD(1.1, 0.1), CD(2.1, 0.0)));
    kv.Add("max_iter", 100);
    kv.Add("eps", 0.001);
    
  }

};
TEST_F(TEST_EvenTemp, Construct) {

  IFactory* factory;
  CreateFactory(kv);
  EXPECT_NO_THROW(factory = CreateFactory(kv));

}
TEST_F(TEST_EvenTemp, BasisSet) {

  IFactory* factory = CreateFactory(kv);

  vector<CSTO>* stos;
  EXPECT_ANY_THROW(factory->GTOSet());
  EXPECT_NO_THROW(stos= factory->STOSet());

  try {
    factory->GTOSet();
  } catch(std::exception& e) {
    cout << e.what() << endl;
  }

  stos = factory->STOSet();

  EXPECT_TRUE(stos != NULL);
  
  EXPECT_EQ(5 + 4, stos->size());
  EXPECT_EQ(1, (*stos)[0].n());
  EXPECT_EQ(2, (*stos)[6].n());
  EXPECT_DOUBLE_EQ(1.2, (*stos)[0].z().real());
  EXPECT_DOUBLE_EQ(1.2 * 1.3, (*stos)[1].z().real());
  EXPECT_DOUBLE_EQ(1.2 * 1.3 * 1.3, (*stos)[2].z().real());
  EXPECT_DOUBLE_EQ(1.1, (*stos)[5].z().real());
  EXPECT_DOUBLE_EQ(1.1*2.1, (*stos)[6].z().real());
  EXPECT_DOUBLE_EQ(1.1*2.1*2.1, (*stos)[7].z().real());
  
}
TEST_F(TEST_EvenTemp, OptTarget) {

  IFactory* factory = CreateFactory(kv);
  HAtomPI<l2func::CSTO>* hatom = NULL;
  hatom = factory->HAtomPiSTO();
  EXPECT_TRUE(hatom != NULL);
  
  CSTO s1(1.1, 2, 1.2);
  CSTO s2(2.1, 2, 1.3);
  EXPECT_NEAR(-0.195858432000038, hatom->OpEle(s1, s2).real(), +0.0000000000001); 
  delete hatom;
  
}
TEST_F(TEST_EvenTemp, Optimizer) {

  IFactory* factory = CreateFactory(kv);
  IOptimizer<CD>* opt(NULL);
  opt = factory->Optimizer();

  EXPECT_TRUE(opt != NULL);
  EXPECT_TRUE(typeid(*opt) == typeid(OptimizerRestricted<CD>));
  EXPECT_EQ(100, opt->max_iter());
  EXPECT_DOUBLE_EQ(0.001, opt->eps());
}
TEST(TEST_OptCBF, optimize) {

  KeysValues kv(",", " ");
  kv.Add<string>("channel", "1s->kp");
  kv.Add<string>("dipole", "length");
  kv.Add<double>("energy", 1.1-0.5);
  kv.Add<string>("basis_type", "STO");
  kv.Add("opt_basis", make_tuple(2, CD(1.4, -0.2)));
  kv.Add("opt_basis", make_tuple(2, CD(0.2, -0.6)));
  kv.Add("max_iter", 100);
  kv.Add("eps", 0.0000001);

  IFactory* factory = CreateFactory(kv);

  IOptTarget* opt_target = factory->OptTarget();
  IOptimizer<CD>*  optimizer = factory->Optimizer();

  VectorXcd zs0(2);
  zs0 << CD(0.8, -0.1), CD(0.4, -0.6);

  OptRes<CD> opt_res = optimizer->Optimize
    (bind(&IOptTarget::Compute, opt_target,
	  _1, _2, _3, _4), zs0);

  EXPECT_TRUE(opt_res.convergence);
  EXPECT_NEAR(0.964095, opt_res.z(0).real(), 0.000004);
  EXPECT_NEAR(-0.0600633, opt_res.z(0).imag(), 0.0000004);
  EXPECT_NEAR(0.664185, opt_res.z(1).real(), 0.000004);
  EXPECT_NEAR(-1.11116, opt_res.z(1).imag(), 0.00004);
}

