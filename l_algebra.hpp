#ifndef L_ALGEBRA_HPP
#define L_ALGEBRA_HPP

#include <Eigen/Core>

using namespace Eigen;

namespace opt_cbf_h {

  // compute aA^(1,0)b
  // aA^(10)b = { \sum(i)  a_kA^(10)_ki b_i   | k }
  template<class F>
  void Calc_a_A10_b(const Matrix<F, Dynamic, 1>& a,
		    const Matrix<F, Dynamic, Dynamic>& A10,
		    const Matrix<F, Dynamic, 1>& b,
		    Matrix<F, Dynamic, 1>* res);
  
  // ab^(1) = { a_k b^(1)_k   | k }
  template<class F>
  void Calc_a_b1(const Matrix<F, Dynamic, 1>& a, 
		 const Matrix<F, Dynamic, 1>& b1, 
		 Matrix<F, Dynamic, 1>* res);

  template<class F>
  void Calc_a_A10_B_A10_b(const Matrix<F, Dynamic, 1>& a, 
			  const Matrix<F, Dynamic, Dynamic>& A10,
			  const Matrix<F, Dynamic, Dynamic>& B,
			  const Matrix<F, Dynamic, 1>& b, 
			  Matrix<F, Dynamic, Dynamic>* res);

  template<class F>
  void Calc_a_A20_b(const Matrix<F, Dynamic, 1>& a, 
		    const Matrix<F, Dynamic, Dynamic>& A20,
		    const Matrix<F, Dynamic, 1>& b, 
		    Matrix<F, Dynamic, Dynamic>* res);

  template<class F>
  void Calc_a_A20_b(const Matrix<F, Dynamic, 1>& a, 
		    const Matrix<F, Dynamic, Dynamic>& A20,
		    const Matrix<F, Dynamic, 1>& b, 
		    Matrix<F, Dynamic, Dynamic>* res);


  

  // compute aA^jb
  // where a and b are vector and A^j is derivative of j th basis:
  // (A^j)_rs = \delta(j,r) (r'|A|s) + \delta(j,s) (r|A|s')
  // input is A10_rs = (r'|A|s)
  template<class F>
  Matrix<F, Dynamic, 1> Calc_a_Aj_b(const Matrix<F, Dynamic, 1>& a,
				    const Matrix<F, Dynamic, Dynamic>& A10,
				    const Matrix<F, Dynamic, 1>& b);

  
  // compute aA^iBA^jb
  // where a and b are vector and A^j is derivative of j th basis:
  // (A^j)_rs = \delta(j,r) (r'|A|s) + \delta(j,s) (r|A|s')
  // input is A10_rs = (r'|A|s)
  template<class F>
  Matrix<F, Dynamic, Dynamic> Calc_a_Ai_B_Aj_b(const Matrix<F, Dynamic, 1>& a,
					       const Matrix<F, Dynamic, Dynamic>& A10,
					       const Matrix<F, Dynamic, Dynamic>& B,
					       const Matrix<F, Dynamic, 1>& b);

  
  // compute (aA^{i,j}a)_ij
  // i == j
  // a_x A^{ii}_xy a_y = a_i A^{20}_iy a_y + a_x A^{02}_xi a_i
  //                   + 2 a_i A^{11}_ii a_i
  template<class F>
  Matrix<F, Dynamic, Dynamic> Calc_a_Aij_a(const Matrix<F, Dynamic, 1>& a,
					   const Matrix<F, Dynamic, Dynamic>& A20,
					   const Matrix<F, Dynamic, Dynamic>& A11);

  // compute (aA^{i,j}b)_xy
  // i == j
  // a_x A^{ii}_xy b_y = a_i A^{20}_iy b_y + a_x A^{02}_xi b_i
  //                   + 2 a_i A^{11}_ii b_i
  template<class F>
  Matrix<F, Dynamic, Dynamic> Calc_a_Aij_b(const Matrix<F, Dynamic, 1>& a,
					   const Matrix<F, Dynamic, Dynamic>& A20,
					   const Matrix<F, Dynamic, Dynamic>& A11,
					   const Matrix<F, Dynamic, 1>& b);


  // ai_A_Bj_b = d_i(a_x) A_xy d_j(B_yz) b_z
  //           = a'_i A_ij B^(1,0)(j, z) b_z +
  //             a'_i A_ix B^(0,1)(x, j) b_j
  template<class F>
  Matrix<F, Dynamic, Dynamic> Calc_ai_A_Bj_b(const Matrix<F, Dynamic, 1>& a1,
					     const Matrix<F, Dynamic, Dynamic>& A,
					     const Matrix<F, Dynamic, Dynamic>& B10,
					     const Matrix<F, Dynamic, 1>& b);


  // compute a^i b = a^i_i * b_i
  template<class F>
  Matrix<F, Dynamic, 1> Calc_ai_b(const Matrix<F, Dynamic, 1>& a1,
				  const Matrix<F, Dynamic, 1>& b);

}

#endif
