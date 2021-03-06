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
#include <timer.hpp>
#include "l_algebra.hpp"
#include "restrict.hpp"
#include "opt.hpp"
#include "opt_cbf.hpp"

using namespace Eigen;
using namespace opt_cbf_h;
using namespace l2func;
using namespace std;
//typedef std::complex<double> CD;

class TwoDimFunc {
  // represent 2d sample function for optimization.
  // f(x, y)    = 0.5 x^4 - 2x^2y + 4y^2 + 8x + 8y
  // dx f(x, y) = 2x^3 -4xy + 8
  // dy f(x, y) = -2x^2 + 8y + 8
  // dxdx f     = 8
  // dydy f     = -4x
  // dxdy f     = -4x
  // the minimum points are:
  // Rule[x, -1.3646556076560374`], Rule[y, -0.5344287681232319`]]
  
public:
  void CalcGradHess(VectorXd z, double* a, VectorXd* g, MatrixXd* h) {
    double x = z(0,0);
    double y = z(1,0);

    *a = 0.5*x*x*x*x - 2.0*x*x*y + 4.0*y*y + 8.0*x + 8.0*y;
    
    (*g)(0, 0) = 2 * x * x * x - 4 * x * y + 8;
    (*g)(1, 0) = -2*x*x + 8*y + 8;

    (*h)(0,0)  = 6 * x * x - 4 * y;
    (*h)(1,0)  = -4 * x;
    (*h)(0,1)  = -4 * x;
    (*h)(1,1)  = 8; 
  }

};
void CalcFuncTest(VectorXd xs, double* a, VectorXd* g,MatrixXd* h) {
  *a = pow(xs(0) - 1.0, 2) + pow(xs(1) -2.0, 2) + 
    pow(xs(2)-xs(3), 2) + pow(xs(3) - 3.0, 2) + pow(xs(4), 2);

  *g = VectorXd::Zero(5);
  (*g) << 
    2 * (xs(0) - 1), 
    2 * (xs(1) - 2),
    2 * (xs(2) - xs(3)),
    -2* (xs(2) - xs(3)) + 2 * (xs(3) -3),
    2 * xs(4);

  *h = MatrixXd::Zero(5, 5);
  *h << 
    2, 0, 0, 0, 0,
    0, 2, 0, 0, 0,
    0, 0, 2,-2, 0,
    0, 0, -2,4, 0,
    0, 0, 0, 0, 2;
}

TEST(first, first) {
  EXPECT_EQ(2, 1 + 1);

  //  typedef Matrix<std::complex<double> > MatrixNcd;
  //  Matrix<std::complex<double> > A =
  MatrixXd mat = MatrixXd::Zero(4, 4);
  MatrixXcd cmat = MatrixXcd::Zero(4,4);
  EXPECT_DOUBLE_EQ(0.0, mat(1,1));
  EXPECT_DOUBLE_EQ(0.0, cmat(1,1).real());

  Matrix<double, Dynamic, 1> v =Matrix<double, Dynamic, 1>::Ones(3);
  EXPECT_DOUBLE_EQ(1.0, v(1));

  Vector3d a; a << 1.0, 2.0, 3.0;
  Vector3d b; b << 0.1, 0.2, 0.3;
  Vector3d c; c = a.array() * b.array();
  EXPECT_DOUBLE_EQ(c(0), 0.1);
}
TEST(first, outer_prod) {

  VectorXd a(3); a << 1.0, 2.0, 3.0;
  VectorXd b(3); b << 0.3, 0.5, 0.6;
  MatrixXd res(3,3); res = a * b.transpose();
  EXPECT_DOUBLE_EQ(0.3, res(0,0));
  
}
TEST(LinearAlgebra, a_Aj_b) {

  // A^0 = (2.0, 2.0, 
  //        2.0, 0.0 )
  // A^1 = (0.0,   1.5)
  //       (1.5,   0.2)
  // aA^0b = (0.1,0.2) (1.6, 0.6) = 0.16 + 0.12 = 0.28
  // aA^1b = (0.1,0.2) (0.75, 0.45+0.1) = 0.075+0.11=0.185

  MatrixXd A10(2,2); A10 << 1.0, 2.0, 1.5, 0.1;
  VectorXd a(2), b(2);  a << 0.1, 0.2;  b << 0.3, 0.5;
  VectorXd res(2);
  Calc_a_Aj_b(a, A10, b, &res);

  EXPECT_DOUBLE_EQ(0.28, res(0));
  EXPECT_DOUBLE_EQ(0.185, res(1));
}
TEST(LinearAlgebra, a_Ai_B_Aj_b) {

   // A^0 = ( 2.0, 2.0,
   //         2.0, 0.0)
   // B   = (1.0, 0.1,
   //        0.1, 0.3)
   // a   = (0.2, 0.4)
   // b   = (0.3, 0.1)

   MatrixXd A10(2,2);   A10 << 1.002, 1.008, 1.003, 1.016;
   VectorXd a(2), b(2); a   << 0.2, 0.4; b << 0.3, 0.1;
   MatrixXd B(2,2);     B   << 1.0, 0.1, 0.1, 0.3;

   MatrixXd res(2,2);
   Calc_a_Ai_B_Aj_b(a, A10, B, b, &res);

   MatrixXd A_d0_d0(2, 2); A_d0_d0 << 2.004, 1.008, 1.008, 0.0;
   EXPECT_DOUBLE_EQ( a.transpose() * A_d0_d0 * B * A_d0_d0 * b,
		     res(0, 0));
 }
TEST(LinearAlgebra, a_Aij_a) {

  MatrixXd A20(2,2); A20 << 6.6, 6.6, 7.2, 7.2;
  MatrixXd A11(2,2); A11 << 1.0, 1.0, 1.0, 1.0;
  VectorXd a(2)    ; a   << 0.2, 0.4;

  MatrixXd res(2,2);
  Calc_a_Aij_a(a, A20, A11, &res);
  EXPECT_DOUBLE_EQ(1.664, res(0,0));
    
}
TEST(LinearAlgebra, real_sto) {

  typedef LinearComb<RSTO> B;
  int num = 3;
  vector<B> us;
  vector<B> d_us;
  vector<B> dd_us;
  for(int i = 0; i < num; i++) {
    double z = 0.1 + i * 0.5;
    RSTO ui = RSTO(2, z, Normalized);
    us.push_back( ui);
    LinearComb<RSTO> d_ui =  D1Normalized(ui);
    LinearComb<RSTO> dd_ui = D2Normalized(ui);
    d_us.push_back(d_ui);
    dd_us.push_back(dd_ui);
  }
  RSTO driv(1.1, 1, 1.0);

  MatrixXd A(num, num);
  MatrixXd A10(num, num);
  MatrixXd A20(num, num);
  MatrixXd A11(num, num);
  MatrixXd A_d0 = MatrixXd::Zero(num, num);
  MatrixXd A_d1 = MatrixXd::Zero(num, num);
  //  MatrixXd A_d00(num, num);
  //  MatrixXd A_d10(num, num);
  VectorXd m(num);
  
  for(int i = 0; i < num; i++) {

    m(i) = CIP(us[i], driv);
    
    for(int j = 0; j < num; j++) {

      A(i, j)   = CIP(us[i],    us[j]);
      A10(i, j) = CIP(d_us[i],  us[j]);
      A20(i, j) = CIP(dd_us[i], us[j]);
      A11(i, j) = CIP(d_us[i], d_us[j]);
    }
  }

  for(int i = 0; i < num; i++) {
    if(i ==0)
      A_d0(0, 0) = 2 * CIP(d_us[0], us[0]);
    else 
      A_d0(0, i) = A_d0(i, 0) = CIP(d_us[0], us[i]);

  
    if(i == 1)
      A_d1(1, 1) = 2 * CIP(d_us[1], us[1]);
    else 
      A_d1(1, i) = A_d1(i, 1) = CIP(d_us[1], us[i]);
  /*
    if(i == 0)
      A_d00(0, 0) = 2 * (CIP(d_us[0], d_us[0]) + CIP(us[0], dd_us[0]));
    else
      A_d00(i, 0) = A_d00(0, i) = 2 * CIP(dd_us[0], us[i]);
    */
  }

  /*
  for(int i = 0; i < num; i++)
    for(int j = 0; j < num; j++) {
      A_d00()
    }
  */

  VectorXd aAia;
  Calc_a_Aj_b(m, A10, m, &aAia);
  double expe =  m.transpose() * A_d0 * m;
  double calc = aAia(0, 0);
  //EXPECT_NEAR(0.0, (calc - expe)/expe, 0.001);
  EXPECT_DOUBLE_EQ(expe, calc);
  EXPECT_NEAR( m.transpose() * A_d1 * m, aAia(1), 0.0001);

  /*
  MatrixXd a_Aij_a;
  Calc_a_Aij_a(m ,A20, A10, &a_Aij_a);
  EXPECT_DOUBLE_EQ( m.transpose() * A_d0 * m, aAia(0, 0));
  */
}
TEST(Optimizer, Newton) {

  IOptimizer<double>* opt = 
    new OptimizerNewton<double>(100, 0.0000001);
  TwoDimFunc func;

  VectorXd x0 = VectorXd::Zero(2);
  x0(0,0) = 1.0; x0(1, 0) = 2.0;
  
  OptRes<double> res = opt->Optimize
    (bind(&TwoDimFunc::CalcGradHess, func, _1, _2, _3, _4), x0);

  EXPECT_EQ(100, opt->max_iter());
  EXPECT_DOUBLE_EQ(0.0000001, opt->eps());
  
  delete opt;
  double eps(0.000000001);

  EXPECT_TRUE(res.convergence);
  EXPECT_TRUE(res.iter_num < 50);
  
  EXPECT_NEAR(res.z(0,0), -1.3646556076560374, eps);
  EXPECT_NEAR(res.z(1,0), -0.5344287681232319, eps);
}
TEST(Optimizer, NewtonWithEventemp) {

  /*
  IRestriction<double>* even_temp = 
    new EvenTemp<double>(2, 1.1, 3.0);
  */
  IRestriction<double>* even_temp = 
    new EvenTemp<double>();
  IOptimizer<double>* opt = 
    new OptimizerRestricted<double>
    (100, 0.000001, even_temp);
  
  VectorXd x0 = VectorXd(2);
  x0 << -1.3, -0.5;
  
  TwoDimFunc func;
  OptRes<double> res = opt->Optimize
    (bind(&TwoDimFunc::CalcGradHess, func, _1, _2, _3, _4), x0);
  
  delete opt;

  double eps(0.000000001);
  EXPECT_TRUE(res.convergence);
  EXPECT_TRUE(res.iter_num < 50);
  EXPECT_NEAR(res.z(0,0), -1.3646556076560374, eps);
  EXPECT_NEAR(res.z(1,0), -0.5344287681232319, eps);  
				     
}
TEST(Optimizer, MultiET) {

  vector<int> idxs(2);
  idxs[0] = 2; idxs[1] = 3;
  IRestriction<double>* et = new MultiEvenTemp<double>(idxs);
  IOptimizer<double>* opt = 
    new OptimizerRestricted<double>(100, 0.00001, et);

  //xs0.resize(5);
    // in this test, the initial guess is choosen as
    // (a,r,b,s) = (0.5, 1.8, 2.2, 0.7)
    // but IRestriction class only support setting by 
    // original sequences.
  VectorXd xs(5);
  xs << 0.5, 0.5*1.8, 2.2, 2.2*0.7, 2.2*0.7*0.7;
  
  OptRes<double> res = opt->Optimize(CalcFuncTest, xs);

  // optimized points are
  // (a,r,b,s) = (1,2,2.48028369537265,0.724491959000516)
  double eps(0.0000001);
  double b0(2.48028369537265);
  double s0(0.724491959000516);
  EXPECT_TRUE(res.convergence);
  EXPECT_NEAR(1.0, res.z(0), eps);
  EXPECT_NEAR(2.0, res.z(1), eps);
  EXPECT_NEAR(b0, res.z(2), eps);
  EXPECT_NEAR(b0*s0, res.z(3), eps);
  EXPECT_NEAR(b0*s0*s0, res.z(4), eps);
}
TEST(Restriction, NoRestriction) {

  VectorXd xs(4);
  xs << 1.1, 2.2, 2.3, 2.4;   
  //IRestriction<double>* no_rest = new NoRestriction<double>(xs);
  IRestriction<double>* no_rest = new NoRestriction<double>();
  no_rest->SetVars(xs);
  EXPECT_EQ(4, no_rest->size());

}
TEST(Restriction, EvenTemp) {
  
  EvenTemp<double> even_temp;
  VectorXd xs(4);
  xs << 1.1, 2.2, 2.3, 2.4; 
  even_temp.SetVars(xs);

  EXPECT_EQ(4, even_temp.size());
  
  EXPECT_DOUBLE_EQ(2.0, even_temp.ratio());
  EXPECT_DOUBLE_EQ(1.1, even_temp.x0());

  xs = even_temp.Xs();
  EXPECT_DOUBLE_EQ(1.1, xs(0));
  EXPECT_DOUBLE_EQ(2.2, xs(1));
  EXPECT_DOUBLE_EQ(4.4, xs(2));

  // consider 2 variable function
  // f(x, y, z)    = x - 1.2y + 2.5sin(z)
  // f(x, rx, rrx) = x - 1.2rx + 2.5sin(rrx)
  // dx f          = 1 - 1.2r +2.5rr cos(rrx)
  // dr f          = -1.2x + 5rx cos(rrx) 

  VectorXd grad_xyz(3);
  grad_xyz << 1.0, -1.2, 2.5 * cos(4.4);
  VectorXd grad = even_temp.Grad(grad_xyz);
  EXPECT_DOUBLE_EQ( 1.0 - 1.2 * 2 + 2.5 * 4 * cos(4.4),
		    grad(0));
  EXPECT_DOUBLE_EQ(-1.2*1.1 + 5 * 2 * 1.1 * cos(4.4),
		   grad(1));
  
  // dxdx f = -2.5rrrr sin(rrx)
  // dxdr f = -1.2 + 5r cos(rrx) - 5rrrx sin(rrx)
  // drdr f = 5x cos(rrx) - 10rrxx sin(rrx)
  MatrixXd hess_xyz(3,3);
  hess_xyz << 
    0.0, 0.0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0, -2.5 * sin(4.4);
  MatrixXd hess = even_temp.Hess(grad_xyz, hess_xyz);
  EXPECT_DOUBLE_EQ(-2.5 * 16 * sin(4.4),
		   hess(0, 0));
  EXPECT_DOUBLE_EQ(-1.2 + 5*2*cos(4.4) - 5*8*1.1*sin(4.4),
		   hess(1,0));
  EXPECT_DOUBLE_EQ(5 * 1.1 * cos(4.4) - 10*4*1.1*1.1*sin(4.4),
		   hess(1,1));
  
}
TEST(Restriction, InRestrictSpace) {
  
  EvenTemp<double> even_temp;
  VectorXd xs0(4);
  xs0 << 1.1, 2.2, 2.3, 2.4; 
  even_temp.SetVars(xs0);

  // now the sequence is 
  // 1.1, 2.2, 4.4, 8.8
  VectorXd xs1 = even_temp.Xs();

  // shift 1.2 for x0 and 3 for ratio
  VectorXd dx(2);
  dx << 1.2, 3;
  even_temp.Shift(dx);

  // now the sequence is
  // 2.3, 11.5, 57.5, ...
  VectorXd xs2 = even_temp.Xs();
  EXPECT_DOUBLE_EQ(2.3, xs2(0));
  EXPECT_DOUBLE_EQ(2.3 * 5.0, xs2(1));
  EXPECT_DOUBLE_EQ(2.3 * 5.0 * 5.0, xs2(2));

}
TEST(Restriction, Interface) {

  IRestriction<double>* even_temp = 
    new EvenTemp<double>();

  VectorXd xs = VectorXd::Zero(3);
  xs << 1.1, 2.2, 2.3;
  even_temp->SetVars(xs);

  EXPECT_EQ(3, even_temp->size());

}
class TEST_MultiET : public ::testing::Test {
public:  

  // this test is based on caluculation in 
  // supply/opt_func.ipynb written in python
  // f(x) : 5 variable function(x=x0,x1,...,x7)
  // f(x) = (x0-1)^2 + (x1-2)^2 + (x2-x3)^2 
  //        (x3-3)^2 x4^2
  // the stationary point for f is 
  // (1,2,3,3,0)
  // 
  // I restricted variable as follows
  // x0 = a, x1 = ar,
  // x2 = b, x3 = bs, x4 = bss
  // the restricted optimization points are
  // (a,r,b,s) = (1,2,2.48028369537265,0.724491959000516)

  // xsTest is location for test (see python notebook)
  VectorXd xsTest;  
  MultiEvenTemp<double>* et;
  virtual void SetUp() {
    
    xsTest.resize(5);
    xsTest << 1,2,2,6,18;

    // idxs represents the number of sequence for each 
    // even-tempered numbers.
    vector<int> idxs(2);
    idxs[0] = 2; idxs[1] = 3;

    et = new MultiEvenTemp<double>(idxs);

  }

};
TEST_F(TEST_MultiET, construct) {
 
  // IRestriction object must set vars before any use.
  et->SetVars(xsTest);
   
  EXPECT_EQ(5, et->size());
  
  EXPECT_EQ(2, get<0>(et->num_x0_r(0)));
  EXPECT_DOUBLE_EQ(1.0, get<1>(et->num_x0_r(0)));
  EXPECT_DOUBLE_EQ(2.0, get<2>(et->num_x0_r(0)));

  EXPECT_EQ(3, get<0>(et->num_x0_r(1)));  
  EXPECT_DOUBLE_EQ(2.0, get<1>(et->num_x0_r(1)));
  EXPECT_DOUBLE_EQ(3.0, get<2>(et->num_x0_r(1)));

  EXPECT_ANY_THROW(et->num_x0_r(-1));
  EXPECT_ANY_THROW(et->num_x0_r(2));
}
TEST_F(TEST_MultiET, Gradient) {
  
  et->SetVars(xsTest);

  // this gradient is calculated in python notebook
  VectorXd grad_normal(5);
  grad_normal << 0, 0, -8, 14, 36;
  VectorXd grad = et->Grad(grad_normal);

  //
  // df_db = 3^0 * (-2) + 3^1*(4) + 3^2*10 
  //       = -2 + 12 + 90 = 100
  EXPECT_DOUBLE_EQ(0.0, grad(0));
  EXPECT_DOUBLE_EQ(0.0, grad(1));
  EXPECT_DOUBLE_EQ(358.0, grad(2));
  EXPECT_DOUBLE_EQ(460.0, grad(3));
  
}
TEST_F(TEST_MultiET, Hess) {

  et->SetVars(xsTest);

  // this hess and grad are calculated in python notebook
  VectorXd grad_normal(5);
  grad_normal << 0, 0, -8, 14, 36;
  MatrixXd hess_normal(5, 5);
  hess_normal << 
    2, 0, 0, 0, 0,
    0, 2, 0, 0, 0,
    0, 0, 2, -2, 0,
    0, 0, -2, 4, 0,
    0, 0, 0, 0, 2;
  MatrixXd hess = et->Hess(grad_normal, hess_normal);

  EXPECT_DOUBLE_EQ(10.0, hess(0,0));
  EXPECT_DOUBLE_EQ(4.0,  hess(1,0));
  EXPECT_DOUBLE_EQ( 0.0, hess(2,0));
  EXPECT_DOUBLE_EQ( 0.0, hess(3,0));

  EXPECT_DOUBLE_EQ(4.0, hess(0,1));
  EXPECT_DOUBLE_EQ(2.0, hess(1,1));
  EXPECT_DOUBLE_EQ(0.0, hess(2,1));
  EXPECT_DOUBLE_EQ(0.0, hess(3,1));

  EXPECT_DOUBLE_EQ(0.0, hess(0,2));
  EXPECT_DOUBLE_EQ(0.0,  hess(1,2));
  EXPECT_DOUBLE_EQ(188.0, hess(2,2));
  EXPECT_DOUBLE_EQ(466.0, hess(3,2));

  EXPECT_DOUBLE_EQ(0.0,   hess(0,3));
  EXPECT_DOUBLE_EQ(0.0,   hess(1,3));
  EXPECT_DOUBLE_EQ(466.0, hess(2,3));
  EXPECT_DOUBLE_EQ(448.0, hess(3,3));

}
TEST_F(TEST_MultiET, Xs) {

  VectorXd xs(5);
  xs << 1,2,2,6,6;

  et->SetVars(xs);
  EXPECT_DOUBLE_EQ (18.0,  et->Xs()(4));

}
TEST_F(TEST_MultiET, Shift) {
  
  // set xs = (1,2,2,6,18)
  et->SetVars(xsTest);
  // (a0, r0, a1, r1) = (1,2,2,3)

  // shift (1, 2, 3, -1);
  VectorXd dx(4);
  dx << 1.0, 2.0, 3.0, -1.0;
  et->Shift(dx);

  // (a0, r0, a1, r1) = (2, 4, 5, 2)
  // xs = (2, 8, 5, 10, 20)
  EXPECT_DOUBLE_EQ(2.0, et->Xs()(0));
  EXPECT_DOUBLE_EQ(8.0, et->Xs()(1));
  EXPECT_DOUBLE_EQ(5.0, et->Xs()(2));
  EXPECT_DOUBLE_EQ(10.0, et->Xs()(3));
  EXPECT_DOUBLE_EQ(20.0, et->Xs()(4));
}
class TestOptSTO: public ::testing::Test {
public:
  VectorXcd zs;
  vector<CSTO> basis_set;
  LinearComb<CSTO> mu_phi;// driven term
  
  IOptTarget* opt_cbf;
  virtual void SetUp() {

    zs = VectorXcd::Zero(2);
    zs(0) = CD(1.4, -0.2);
    zs(1) = CD(0.2, -0.6);

    basis_set.resize(2);
    for(int i = 0; i < 2; i++)
      basis_set[i] = CSTO(2, zs(i), Normalized);

    mu_phi += 1.0 * CSTO(CD(2.0, 0.0), 2, CD(1.0, 0.0));

    // w = 1.1, E0 = -0.5  =>  E = 0.6
    HLikeAtom<CD> h_atom(2, 1.0, 1);
    opt_cbf = new OptCBF<CSTO, CSTO>(basis_set, h_atom, mu_phi, 0.6);
  }
  virtual ~TestOptSTO() {
    //delete h_atom;
    //delete opt_cbf;
  }
};
TEST_F(TestOptSTO, matrix) {

  /*
  MatrixXcd D00, D10, D20, D11;
  opt_cbf->computeMatrix(&D00, &D10, &D20, &D11);

  double eps(0.000000001);
  EXPECT_NEAR(-0.34, D00(0,0).real(), eps);
  EXPECT_NEAR(-0.18, D00(0,0).imag(), eps);
  EXPECT_NEAR(+0.777278, D00(1,0).real(),
	      +0.000001);
  EXPECT_NEAR(-0.987452, D00(1,0).imag(),
	      +0.000001);
  EXPECT_NEAR(-0.86, D00(1,1).real(), eps);
  EXPECT_NEAR(+0.18, D00(1,1).imag(), eps);  
  */
  
}
TEST_F(TestOptSTO, vector) {
  
  /*
  VectorXcd m0, m1, m2;
  opt_cbf->computeVector(&m0, &m1, &m2);
  EXPECT_NEAR(1.62413 ,m0(0).real(), 0.00001);
  EXPECT_NEAR(0.0991355 ,m0(0).imag(), 0.0000001);
  EXPECT_NEAR(-2.81312, m0(1).real(), 0.00001);
  EXPECT_NEAR(2.92198 , m0(1).imag(), 0.00001);

  EXPECT_NEAR(-0.525733, m1(0).real(), 0.000001);
  EXPECT_NEAR(0.0943892, m1(0).imag(), 0.000001);
  EXPECT_NEAR(-0.226781, m1(1).real(), 0.000001);
  EXPECT_NEAR(-11.9481, m1(1).imag(), 0.0001);  
  */

}
TEST_F(TestOptSTO, Construct) {
  
  EXPECT_DOUBLE_EQ(1.4, zs[0].real());
  
}
TEST_F(TestOptSTO, alpha_grad_hess) {

  CD alpha;
  VectorXcd grad;
  MatrixXcd hess;

  opt_cbf->Compute(zs, &alpha, &grad, &hess);

  double eps(0.000000001);
  EXPECT_NEAR(-6.930263229774504, alpha.real(), eps);
  EXPECT_NEAR(0.9284529898976182, alpha.imag(), eps);

  EXPECT_NEAR(-8.61257434411908, grad(0, 0).real(),    0.00001);
  EXPECT_NEAR(+0.02011898864080297, grad(0, 0).imag(), 0.00001);
  EXPECT_NEAR(+4.725266440501804, grad(1, 0).real(),   0.00001);
  EXPECT_NEAR(0.7486897022445742, grad(1, 0).imag(),   0.00001);

  eps = 0.00001;
  EXPECT_NEAR(-22.6814, hess(0, 0).real(),
	        0.0001);
  EXPECT_NEAR(-2.49483, hess(0, 0).imag(),
	      +0.00001);
  EXPECT_NEAR(27.1967, hess(1, 0).real(),
	      +0.0001);
  EXPECT_NEAR(9.97019, hess(1, 0).imag(),
	      0.00001);
  EXPECT_NEAR(27.1967, hess(0, 1).real(),
	      +0.0001);
  EXPECT_NEAR(9.97019, hess(0, 1).imag(),
	      0.00001);
  EXPECT_NEAR(-11.6649, hess(1, 1).real(),
	      +0.0001);
  EXPECT_NEAR(-2.89574, hess(1, 1).imag(),
	      +0.00001);
}
TEST_F(TestOptSTO, optimization) {

  VectorXcd zs0(2);
  zs0 << CD(0.8, -0.1), CD(0.4, -0.6);

  IOptimizer<CD>* opt = new OptimizerNewton<CD>(100, 0.00001);
  OptRes<CD> opt_res = opt->Optimize
    (bind(&IOptTarget::Compute, opt_cbf,
	  _1, _2, _3, _4), zs0);

  //  opt_target->WritePsi("supply/psi.out", 10.0, 1.0);

//  IOptimizer<CD>* opt = new OptimizerNewton<CD>();
//  OptRes<complex<double> > opt_res = opt->Optimize
//    (bind(&OptCBF<CSTO>::Compute, opt_cbf, _1, _2, _3, _4),zs0);

  EXPECT_TRUE(opt_res.convergence);
  EXPECT_NEAR(0.964095, opt_res.z(0).real(), 0.000001);
  EXPECT_NEAR(-0.0600633, opt_res.z(0).imag(), 0.0000001);
  EXPECT_NEAR(0.664185, opt_res.z(1).real(), 0.000001);
  EXPECT_NEAR(-1.11116, opt_res.z(1).imag(), 0.00001);
}
TEST(TestOptCbf, TimeCheck) {

  // time comparing
  // non-opt : 0.261512 s
  // symmetry: 0.204394 s
  // sym + L : 0.069230 s
  
  int num(100);
  VectorXcd zs = VectorXcd::Zero(num);
  for(int i = 0; i < num; i++) 
    zs(i) = CD(0.1 * i, -0.05 * i);
  
  vector<CSTO> basis_set(num);
  for(int i = 0; i < num; i++)
    basis_set[i] = CSTO(2, zs(i), Normalized);

  LinearComb<CSTO> mu_phi;
  mu_phi += 1.0 * CSTO(CD(2.0, 0.0), 2, CD(1.0, 0.0));
  HLikeAtom<CD> h_atom(2, 1.0, 1);

  OptCBF<CSTO,CSTO>* opt_cbf = 
    new OptCBF<CSTO, CSTO>(basis_set, h_atom, mu_phi, 0.6);

  CD alpha;
  VectorXcd grad;
  MatrixXcd hess;
  Timer timer;
  
  timer.Start("calc");
  opt_cbf->Compute(zs, &alpha, &grad, &hess);
  timer.End("calc");

  cout << "Time: " << timer.GetTime("calc") << endl;
  

}



