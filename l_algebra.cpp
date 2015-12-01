#include "l_algebra.hpp"

namespace opt_cbf_h {

  template<class F>
  void Calc_a_A10_b(const Matrix<F, Dynamic, 1>& a,
		    const Matrix<F, Dynamic, Dynamic>& A10,
		    const Matrix<F, Dynamic, 1>& b,
		    Matrix<F, Dynamic, 1>* res) {

    Matrix<F, Dynamic, 1> A10_b = A10 * b;
    *res = a.array() * A10_b.array() ;

  }
  template void Calc_a_A10_b(const VectorXcd&, const MatrixXcd&, const VectorXcd&,
			     VectorXcd*);  

  template<class F>
  void Calc_a_b1(const Matrix<F, Dynamic, 1>& a, 
		 const Matrix<F, Dynamic, 1>& b1, 
		 Matrix<F, Dynamic, 1>* res) {
    *res = a.array() * b1.array();
  }
  template void Calc_a_b1(const VectorXcd&, const VectorXcd&, VectorXcd*);

  


  template<class F>
  Matrix<F, Dynamic, 1> Calc_a_Aj_b(const Matrix<F, Dynamic, 1>& a,
				    const Matrix<F, Dynamic, Dynamic>& A10,
				    const Matrix<F, Dynamic, 1>& b) {

    Matrix<F, Dynamic, 1> A10_b = A10 * b;
    Matrix<F, Dynamic, 1> res1  = a.array() * A10_b.array() ;

    Matrix<F, Dynamic, 1> A10_a = A10 * a;
    Matrix<F, Dynamic, 1> res2 = b.array() * A10_a.array();

    return res1 + res2;
  }

  template<class F>
  Matrix<F, Dynamic, Dynamic> Calc_a_Ai_B_Aj_b(const Matrix<F, Dynamic, 1>& a,
					       const Matrix<F, Dynamic, Dynamic>& A10,
					       const Matrix<F, Dynamic, Dynamic>& B,
					       const Matrix<F, Dynamic, 1>& b) {

    typedef Matrix<F, Dynamic, Dynamic> M;

    M tmp1 = (a * (A10*b).transpose());
    M tmp2 = A10 * B;
    M t1 = tmp1.array() * tmp2.array();

    M t2 = ((A10*a) * (A10*b).transpose()).array() * B.array();
    M t3 = (A10*B*A10.transpose()).array() * (a*b.transpose()).array();
    M t4 = (A10*a * b.transpose() ).array()
	    * (B*A10.transpose()).array();
    return t1 + t2 + t3 + t4;

  }

  template<class F>
  Matrix<F, Dynamic, Dynamic> Calc_a_Aij_a(const Matrix<F, Dynamic, 1>& a,
					   const Matrix<F, Dynamic, Dynamic>& A20,
					   const Matrix<F, Dynamic, Dynamic>& A11){
    Matrix<F, Dynamic, Dynamic> res;
    res = F(2) * ((a * a.transpose()).array() * A11.array());
    res+=F(2) * ((A20 * a).array() * a.array()).matrix().asDiagonal();
    return res;
  }

  template<class F>
  Matrix<F, Dynamic, Dynamic> Calc_a_Aij_b(const Matrix<F, Dynamic, 1>& a,
					   const Matrix<F, Dynamic, Dynamic>& A20,
					   const Matrix<F, Dynamic, Dynamic>& A11,
					   const Matrix<F, Dynamic, 1>& b) {
    Matrix<F, Dynamic, Dynamic> res;
    res = (a * b.transpose()).array() * A11.array();
    res+= ((b * a.transpose()).array() * A11.array()).matrix();
    res+= ((A20 * a).array() * b.array()).matrix().asDiagonal();
    res+= ((A20 * b).array() * a.array()).matrix().asDiagonal();
    return res;

  }


  template<class F>
  Matrix<F, Dynamic, Dynamic> Calc_ai_A_Bj_b(const Matrix<F, Dynamic, 1>& a1,
					     const Matrix<F, Dynamic, Dynamic>& A,
					     const Matrix<F, Dynamic, Dynamic>& B10,
					     const Matrix<F, Dynamic, 1>& b) {

    return 
      (a1 * (B10 * b).transpose()).array() * A.array() +
      (a1 * b.transpose()).array() * (A * B10.transpose()).array();
  }

  template<class F>
  Matrix<F, Dynamic, 1> Calc_ai_b(const Matrix<F, Dynamic, 1>& a1,
				  const Matrix<F, Dynamic, 1>& b) {

    return  a1.array() * b.array();
  }

  typedef std::complex<double> CD;
  template VectorXcd Calc_a_Aj_b<CD>(const VectorXcd&,
				     const MatrixXcd&,
				     const VectorXcd&);
  template MatrixXcd Calc_a_Ai_B_Aj_b<CD>(const VectorXcd&,
					  const MatrixXcd&,
					  const MatrixXcd&,
					  const VectorXcd&);
  template MatrixXcd Calc_a_Aij_a<CD>(const VectorXcd&,
				      const MatrixXcd&,
				      const MatrixXcd&);
  template MatrixXcd Calc_a_Aij_b<CD>(const VectorXcd& ,
				      const MatrixXcd&,
				      const MatrixXcd&,
				      const VectorXcd&);
  template MatrixXcd Calc_ai_A_Bj_b<CD>(const VectorXcd&,
					const MatrixXcd&,
					const MatrixXcd&,
					const VectorXcd&);
  template VectorXcd Calc_ai_b<CD>(const VectorXcd&,
				   const VectorXcd&);

}
