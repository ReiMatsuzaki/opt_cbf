#include <iostream>
#include <complex>
#include <string>
#include <vector>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <gtest/gtest.h>
#include <Eigen/Core>
#include <l2func.hpp>
#include <keys_values.hpp>
#include "l_algebra.hpp"
#include "restrict.hpp"
#include "opt.hpp"
#include "opt_cbf.hpp"
#include "factory.hpp"

using namespace Eigen;
using namespace opt_cbf_h;
using namespace l2func;
using namespace std;

OptRes<CD> RunOptimize(const KeysValues& kv) {

  IFactory* factory = CreateFactory(kv);
  IOptTarget* opt_target = factory->OptTarget();
  IOptimizer<CD>*  optimizer = factory->Optimizer();

  int num = factory->BasisSize();
  VectorXcd zs0(num);
  factory->GetZs(&zs0);

  OptRes<CD> opt_res = optimizer->Optimize
    (bind(&IOptTarget::Compute, opt_target,
	  _1, _2, _3, _4), zs0);  

  return opt_res;
}

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
  //EXPECT_ANY_THROW(dynamic_cast<FactoryMono<CSTO,CSTO>* > (factory));
  FactoryMono<CSTO,CSTO>* ptr =
    dynamic_cast<FactoryMono<CSTO,CSTO>* > (factory);
  EXPECT_TRUE(ptr == NULL);

  gtos = (dynamic_cast<FactoryMono<CGTO,CSTO>* >(factory))->BasisSet();

  EXPECT_TRUE(gtos != NULL);
  
  EXPECT_EQ(3, gtos->size());
  EXPECT_EQ(2, (*gtos)[1].n());
  EXPECT_DOUBLE_EQ(1.4, (*gtos)[2].z().real());    
    
}
TEST_F(TEST_MonoBasis, Getter) {

  IFactory* factory = CreateFactory(kv);
  VectorXcd zs(3);
  factory->GetZs(&zs);
  EXPECT_DOUBLE_EQ(1.2, zs(0).real());
  EXPECT_DOUBLE_EQ(3.4, zs(1).real());
  EXPECT_DOUBLE_EQ(1.4, zs(2).real());

  EXPECT_EQ(3, factory->BasisSize());

}
TEST_F(TEST_MonoBasis, OptTarget) {

  IFactory* factory = CreateFactory(kv);
  EXPECT_TRUE(factory != NULL);
  IOptTarget* target = factory->OptTarget();
  EXPECT_TRUE(target != NULL);
  //HAtomPI<l2func::CGTO>* hatom = NULL;
  
  //  hatom = factory->HAtomPiGTO();

  //  EXPECT_TRUE(hatom != NULL);

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

  vector<CSTO>* stos = 
    dynamic_cast<FactoryEvenTemp<CSTO,CSTO>*>(factory)->BasisSet();
  
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
TEST_F(TEST_EvenTemp, Getter) {

  IFactory* factory = CreateFactory(kv);
  VectorXcd zs(9);

  CD za(1.2, 0.1), ra(1.3, 0.0), zb(1.1, 0.1), rb(2.1, 0.0);

  factory->GetZs(&zs);
  EXPECT_DOUBLE_EQ((za).real(),    zs(0).real());
  EXPECT_DOUBLE_EQ((za*ra).real(), zs(1).real());
  EXPECT_DOUBLE_EQ((za*ra*ra).real(), zs(2).real());
  EXPECT_DOUBLE_EQ((zb*rb).real(), zs(6).real());

  EXPECT_EQ(9, factory->BasisSize());

}
TEST_F(TEST_EvenTemp, OptTarget) {

  IFactory* factory = CreateFactory(kv);
  IOptTarget* target(NULL);
  target = factory->OptTarget();
  EXPECT_TRUE(target != NULL);

  /*  
  CSTO s1(1.1, 2, 1.2);
  CSTO s2(2.1, 2, 1.3);

  EXPECT_NEAR(-0.195858432000038, hatom->OpEle(s1, s2).real(), +0.0000000000001); 
  delete hatom;
  */
  
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
TEST(TEST_OptCBF, 2sto) {

  KeysValues kv(",", " ");
  kv.Add<string>("channel", "1s->kp");
  kv.Add<string>("dipole", "length");
  kv.Add<double>("energy", 1.1-0.5);
  kv.Add<string>("basis_type", "STO");
  kv.Add("opt_basis", make_tuple(2, CD(0.8, -0.1)));
  kv.Add("opt_basis", make_tuple(2, CD(0.4, -0.6)));
  kv.Add("max_iter", 100);
  kv.Add("eps", 0.0000001);

  OptRes<CD> opt_res = RunOptimize(kv);

  EXPECT_TRUE(opt_res.convergence);
  EXPECT_NEAR(0.964095, opt_res.z(0).real(), 0.000004);
  EXPECT_NEAR(-0.0600633, opt_res.z(0).imag(), 0.0000004);
  EXPECT_NEAR(0.664185, opt_res.z(1).real(), 0.000004);
  EXPECT_NEAR(-1.11116, opt_res.z(1).imag(), 0.00004);
}
TEST(TEST_OptCBF, et_2sto) {

  KeysValues kv(",", " ");
  kv.Add<string>("channel", "1s->kp");
  kv.Add<string>("dipole", "length");
  kv.Add<double>("energy", 1.1-0.5);
  kv.Add<string>("basis_type", "STO");
  kv.Add<int>("max_iter", 100);
  kv.Add<double>("eps", 0.00001);

  kv.Add("opt_et_basis", make_tuple(2, 2, CD(0.96, -0.06), CD(0.6,-1.0)));

  OptRes<CD> opt_res = RunOptimize(kv);
  EXPECT_TRUE(opt_res.convergence);
  EXPECT_NEAR(0.964095, opt_res.z(0).real(), 0.000004);
  EXPECT_NEAR(-0.0600633, opt_res.z(0).imag(), 0.0000004);
  EXPECT_NEAR(0.664185, opt_res.z(1).real(), 0.000004);
  EXPECT_NEAR(-1.11116, opt_res.z(1).imag(), 0.00004);

  /*
convergence: true
iter_num: 4
opt_et_basis: 2 2 (0.964095,-0.0600633) 
coef: (-2.00517,0.600012)
coef: (0.425609,-0.0661265)
alpha: (-4.86381,0.755562)
read_time: 0.000519
calc_time: 0.0027
  */
  
}
TEST(TEST_OptCBF, et_sto_et_sto) {

  KeysValues kv(",", " ");
  kv.Add<string>("channel", "1s->kp");
  kv.Add<string>("dipole", "length");
  kv.Add<double>("energy", 0.5);
  kv.Add<string>("basis_type", "STO");
  kv.Add<int>("max_iter", 300);
  kv.Add<double>("eps", 0.000001);

  kv.Add("opt_et_basis", make_tuple(2, 2, CD(0.9928, -0.0037), CD(0.996,-0.416)));
  kv.Add("opt_et_basis", make_tuple(2, 2, CD(0.6708, -0.9126), CD(0.879,-0.418)));

  // opt zs:
  //complex(0.992828,-0.00370088),
  //complex(0.987316,-0.417238),
  //complex(0.670764,-0.912606),
  //complex(0.207802,-1.08305)]
  //
  // its ratio
  // 0.9960008804285743-0.41653934041106766j),
  // (0.9078663907455836-0.5406674690353385j),
  // (0.8791758959357688-0.41849115682064325j)]

  double eps(0.0001);
  OptRes<CD> opt_res = RunOptimize(kv);
  EXPECT_TRUE(opt_res.convergence);
  EXPECT_NEAR(0.992828, opt_res.z(0).real(),eps);
  EXPECT_NEAR(-0.00370088, opt_res.z(0).imag(), eps);
  EXPECT_NEAR(0.987316, opt_res.z(1).real(), eps);
  EXPECT_NEAR(-0.417238, opt_res.z(1).imag(), eps);
  EXPECT_NEAR(0.670764, opt_res.z(2).real(), eps);
  EXPECT_NEAR(-0.912606, opt_res.z(2).imag(), eps);
  EXPECT_NEAR(0.207802, opt_res.z(3).real(), eps);
  EXPECT_NEAR(-1.08305, opt_res.z(3).imag(), eps);

}

