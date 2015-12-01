#include <gtest/gtest.h>
#include <l2func.hpp>
#include "l_algebra.hpp"
#include "add_gtest.hpp"
#include "opt_target.hpp"
#include "opt_target_impl.hpp"

using namespace Eigen;
using namespace opt_cbf_h;
using namespace l2func;

TEST(LinearAlgebra, a_Aj_b) {

  // A^0 = (2.0, 2.0, 
  //        2.0, 0.0 )
  // A^1 = (0.0,   1.5)
  //       (1.5,   0.2)
  // aA^0b = (0.1,0.2) (1.6, 0.6) = 0.16 + 0.12 = 0.28
  // aA^1b = (0.1,0.2) (0.75, 0.45+0.1) = 0.075+0.11=0.185

  MatrixXcd A10(2,2); A10 << 1.0, 2.0, 1.5, 0.1;
  VectorXcd a(2), b(2);  a << 0.1, 0.2;  b << 0.3, 0.5;
  VectorXcd res(2);
  res = Calc_a_Aj_b(a, A10, b);

  EXPECT_C_EQ(0.28, res(0));
  EXPECT_C_EQ(0.185, res(1));
}
TEST(LinearAlgebra, a_A01_b) {

  MatrixXcd A10(2,2); A10 << 1.0, 2.0, 1.5, 0.1;
  VectorXcd a(2), b(2);  a << 0.1, 0.2;  b << 0.3, 0.5;

  VectorXcd res(2);
  VectorXcd res0(2), res1(2);

  res = Calc_a_Aj_b(a, A10, b);
  Calc_a_A10_b(a, A10, b, &res0);
  Calc_a_A10_b(b, A10, a, &res1);

  EXPECT_C_EQ(res(0), res1(0) + res0(0));
  EXPECT_C_EQ(res(1), res1(1) + res0(1));
}
TEST(LinearAlgebra, a_b1) {

  VectorXcd a(2); a << 1.0, 2.0;
  VectorXcd b1(2); b1 << 3.0, -2.0;
  VectorXcd res(2);
  Calc_a_b1(a, b1, &res);
  EXPECT_C_EQ(3.0, res(0));
  EXPECT_C_EQ(-4.0, res(1));

}
TEST(LinearAlgebra, sym) {

  MatrixXcd A10(2,2); A10 << 1.0, 2.0, 1.5, 0.1;
  MatrixXcd A20(2,2); A20 << 1.0, 2.0, 1.5, 0.1;
  MatrixXcd A11(2,2); A11 << 1.0, 2.0, 2.0, 0.1;
  VectorXcd a(2), b(2);  a << 0.1, 0.2;  b << 0.3, 0.5;

  MatrixXcd res = Calc_a_Aij_b(a, A20, A11, b);

  EXPECT_C_EQ(0.0, (res - res.transpose())(0, 1));

  res = Calc_a_Aij_b(a, A20, A11, b);
  EXPECT_C_EQ(0.0, (res - res.transpose())(0, 1));

  res = Calc_a_Ai_B_Aj_b(a, A10, A11, b) + Calc_a_Ai_B_Aj_b(b, A10, A11, a);
  EXPECT_C_EQ(0.0, (res - res.transpose())(0, 1));
  

}
class TestOptSTO: public ::testing::Test {
public:
  IOptTarget* opt_cbf;
  virtual void SetUp() {

    HLikeAtom<CD> hatom;
    HLength<1, 0, 1, CD> mu_phi(hatom);
    HminusEOp<1,CD> lop(hatom, 0.6);
    
    std::vector<NormalCSTO> basis_set;
    basis_set.push_back(NormalCSTO(2, CD(1.4, -0.2)));
    basis_set.push_back(NormalCSTO(2, CD(0.2, -0.6)));
    
    OptSTOLength *opt_target = new OptSTOLength(mu_phi.value, 
						basis_set, 
						mu_phi.value, 
						lop.value);
    opt_cbf = opt_target;
  }
  virtual void TearDown() {
    delete opt_cbf;
  }
};

TEST_F(TestOptSTO, AlphaGrad) {

  CD alpha;
  VectorXcd grad;
  MatrixXcd hess;

  VectorXcd zs(2); zs << CD(1.4, -0.2), CD(0.2, -0.6);

  opt_cbf->Compute(zs, &alpha, &grad, &hess);
  double eps;
  
  eps = 0.00000001;
  EXPECT_C_NEAR(CD(-6.930263229774504, 0.9284529898976182), alpha, eps);

  eps = 0.00001;
  EXPECT_C_NEAR(CD(-8.61257434411908,  0.02011898864080297), grad(0, 0), eps);
  EXPECT_C_NEAR(CD(+4.725266440501804, 0.7486897022445742),  grad(1, 0), eps);

  eps = 0.0001;
  EXPECT_C_NEAR(CD(-22.6814, -2.49483), hess(0, 0), eps);
  EXPECT_C_NEAR(CD(27.1967, 9.97019),   hess(1, 0), eps);
  EXPECT_C_NEAR(hess(0, 1),   hess(1, 0), pow(10.0, -10.0));
  EXPECT_C_NEAR(CD(-11.6649, -2.89574), hess(1, 1), eps);
  
}

int main (int argc, char **args) {
  ::testing::InitGoogleTest(&argc, args);
  return RUN_ALL_TESTS();
}
