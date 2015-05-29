#include <iostream>
#include <complex>
#include <string>
#include <vector>
#include <tr1/functional>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <l2func.hpp>
#include "l_algebra.hpp"
#include "opt.hpp"
#include "driv.hpp"
#include "opt_cbf.hpp"

using namespace Eigen;
using namespace opt_cbf_h;
using namespace l2func;
using namespace std::tr1;
using namespace std::tr1::placeholders;
using namespace std;
typedef std::complex<double> CD;

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
  void CalcGradHess(VectorXd z, VectorXd* g, MatrixXd* h) {
    double x = z(0,0);
    double y = z(1,0);
    
    (*g)(0, 0) = 2 * x * x * x - 4 * x * y + 8;
    (*g)(1, 0) = -2*x*x + 8*y + 8;

    (*h)(0,0)  = 6 * x * x - 4 * y;
    (*h)(1,0)  = -4 * x;
    (*h)(0,1)  = -4 * x;
    (*h)(1,1)  = 8; 
  }
};

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

  int num = 3;
  vector<RSTO> us(num);
  vector<RSTO> d_us(num);
  vector<RSTO> dd_us(num);
  for(int i = 0; i < num; i++) {
    double z = 0.1 + i * 0.5;
    us[i]   = RSTO(          1.0, 2, z);
    d_us[i] = DBasis<1,RSTO>(1.0, 2, z);
    dd_us[i]= DBasis<2,RSTO>(1.0, 2, z);
  }
  RSTO driv(1.1, 1, 1.0);

  MatrixXd A(num, num);
  MatrixXd A10(num, num);
  MatrixXd A20(num, num);
  MatrixXd A11(num, num);
  MatrixXd A_d0(num, num);
  MatrixXd A_d1(num, num);
  MatrixXd A_d00(num, num);
  MatrixXd A_d10(num, num);
  VectorXd m(num);
  
  for(int i = 0; i < num; i++) {

    m(i, 0) = CIP(us[i], driv);
    
    for(int j = 0; j < num; j++) {

      A(i, j)   = CIP(us[i],    us[j]);
      A10(i, j) = CIP(d_us[i],  us[j]);
      A20(i, j) = CIP(d_us[i], d_us[j]);
      A11(i, j) = CIP(d_us[i], d_us[j]);
    }
  }

  for(int i = 0; i < num; i++) {
    if(i ==0)
      A_d0(i, 0) = 2 * CIP(d_us[0], us[0]);
    else 
      A_d0(0, i) = A_d0(i, 0) = CIP(d_us[0], us[i]);

    if(i == 1)
      A_d1(1, 1) = 2 * CIP(d_us[1], us[1]);
    else 
      A_d1(1, i) = A_d1(i, 1) = CIP(d_us[i], us[1]);

    if(i == 0)
      A_d00(0, 0) = 2 * (CIP(d_us[0], d_us[0]) + CIP(us[0], dd_us[0]));
    else
      A_d00(i, 0) = A_d00(0, i) = 2 * CIP(dd_us[0], us[i]);
  }

  /*
  for(int i = 0; i < num; i++)
    for(int j = 0; j < num; j++) {
      A_d00()
    }
  */

  VectorXd aAia;
  Calc_a_Aj_b(m, A10, m, &aAia);
  EXPECT_DOUBLE_EQ( m.transpose() * A_d0 * m, aAia(0, 0));
  EXPECT_DOUBLE_EQ( m.transpose() * A_d1 * m, aAia(1, 0));

  MatrixXd a_Aij_a;
  Calc_a_Aij_a(m ,A20, A10, &a_Aij_a);
  EXPECT_DOUBLE_EQ( m.transpose() * A_d0 * m, aAia(0, 0));
		    
}
TEST(Optimizer, Newton) {

  IOptimizer<double>* opt = new OptimizerNewton<double>();
  TwoDimFunc func;


  VectorXd x0 = VectorXd::Zero(2);
  x0(0,0) = 1.0; x0(1, 0) = 2.0;
  
  OptRes<double> res = opt->Optimize
    (bind(&TwoDimFunc::CalcGradHess, func, _1, _2, _3), x0);
  
  delete opt;
  double eps(0.000000001);

  EXPECT_TRUE(res.convergence);
  EXPECT_TRUE(res.iter_num < 50);
  
  EXPECT_NEAR(res.z(0,0), -1.3646556076560374, eps);
  EXPECT_NEAR(res.z(1,0), -0.5344287681232319, eps);
}
TEST(Driv, Construct) {

  LinearComb<RSTO> mu_phi;
  mu_phi += 1.0 * RSTO(1.0, 2, 1.0);

  IDrivSystem<RSTO>* hatom = new HAtomPI<RSTO>(1, 1.0, 0.5, mu_phi);
  // HAtomPI<double>* hatom = new HAtomPI<double>(1, 1.0, 0.5, mu_phi);
  RSTO s1(1.1, 2, 1.2);
  RSTO s2(2.1, 2, 1.3);
  double h01 = hatom->OpEle(s1, s2);
  EXPECT_NEAR(-0.195858432000038, h01,
	      +0.0000000000001);
}
TEST(STL, vector) {

  vector<int> xs(10, 42);
  vector<int>::const_iterator it0 = xs.begin() + 1;
  vector<int>::const_iterator it1 = xs.begin() + 3;
  vector<int> ys(it0, it1); // ys have copied value

  ys[0] = 44;
  ys[1] = 45;

  EXPECT_EQ(xs[0], 42);
  EXPECT_EQ(xs[1], 42);
  EXPECT_EQ(xs[2], 42);
  EXPECT_EQ(xs[3], 42);
  EXPECT_EQ(xs[4], 42);
  EXPECT_EQ(2, ys.size());
  EXPECT_EQ(44, ys[0]);

  EXPECT_EQ(2, std::distance(it0, it1));
}
TEST(OptCBF, construct) {

  vector<CSTO> basis_set(2);
  basis_set[0] = CSTO(2, CD(2.5, -0.4), Normalized);
  basis_set[1] = CSTO(2, CD(1.5, -0.9), Normalized);
  
  LinearComb<CSTO> mu_phi;
  mu_phi += 1.0 * RSTO(1.0, 2, 1.0);
  
  HAtomPI<CSTO> h_atom(1, 1.0, 1.1, mu_phi);
  
  OptCBF<CSTO> opt_cbf(basis_set, h_atom);

  vector<std::complex<double> > zs(2);
  zs[0] = 2.1; zs[1] = 2.2;
  VectorXcd grad;
  MatrixXcd hess;
  std::complex<double> alpha;
  opt_cbf.Compute(zs, &alpha, &grad, &hess);
  cout << alpha << endl;
  
}

