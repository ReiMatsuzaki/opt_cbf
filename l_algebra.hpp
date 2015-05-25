#ifndef L_ALGEBRA_HPP
#define L_ALGEBRA_HPP

#include <Eigen/Core>

using namespace Eigen;

namespace opt_cbf_h {

  // compute aA^jb
  // where a and b are vector and A^j is derivative of j th basis:
  // (A^j)_rs = \delta(j,r) (r'|A|s) + \delta(j,s) (r|A|s')
  // input is A10_rs = (r'|A|s)
  template<class F>
  void Calc_a_Aj_b(const Matrix<F, Dynamic, 1>& a,
		   const Matrix<F, Dynamic, Dynamic>& A10,
		   const Matrix<F, Dynamic, 1>& b,
		   Matrix<F, Dynamic, 1>* res) {

    Matrix<F, Dynamic, 1> A10_b = A10 * b;
    Matrix<F, Dynamic, 1> res1  = a.array() * A10_b.array() ;

    Matrix<F, Dynamic, 1> A10_a = A10 * a;
    Matrix<F, Dynamic, 1> res2 = b.array() * A10_a.array();

    *res = res1 + res2;
  }
  
  template<class F>
  void Calc_a_Ai_B_Aj_b(const Matrix<F, Dynamic, 1>& a,
			const Matrix<F, Dynamic, Dynamic>& A10,
			const Matrix<F, Dynamic, Dynamic>& B,
			const Matrix<F, Dynamic, 1>& b,
			Matrix<F, Dynamic, Dynamic>* res) {

    typedef Matrix<F, Dynamic, Dynamic> M;

    M tmp1 = (a * (A10*b).transpose());
    M tmp2 = A10 * B;
    M t1 = tmp1.array() * tmp2.array();

    

    M t2 = ((A10*a) * (A10*b).transpose()).array() * B.array();
    M t3 = (A10*B*A10.transpose()).array() * (a*b.transpose()).array();
    M t4 = (A10*a * b.transpose() ).array()
	    * (B*A10.transpose()).array();
    *res = t1 + t2 + t3 + t4;
    /*
      ((A10 * a) * (A10*b).transpose()).array() * B.array() +
      (A10*B*A10.transpose()).transpose()).array() * (a*b.transpose()).array() +
      ((A10*a)*b.transpose()).array() * (B*A10.transpose()).array();
    */
    //    *res = t1 + t2 + t3 + t4;
  }

  template<class F>
  void Calc_a_Aij_a(const Matrix<F, Dynamic, 1>& a,
		    const Matrix<F, Dynamic, Dynamic>& A20,
		    const Matrix<F, Dynamic, Dynamic>& A11,
		    Matrix<F, Dynamic, Dynamic>* res) {
    *res =
      F(2) * ((a * a.transpose()).array() * A11.array());

    int n_row = res->rows();
    Matrix<F, Dynamic, 1> vec = F(2) * a.array() * (A20*a).array();
    
    for(int i = 0; i < n_row; i++) 
      (*res)(i) += vec(i);
  }
}

#endif
