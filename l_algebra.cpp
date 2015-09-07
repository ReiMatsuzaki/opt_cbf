#include "l_algebra.hpp"

namespace opt_cbf_h {
  
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

  }

  template<class F>
  void Calc_a_Aij_a(const Matrix<F, Dynamic, 1>& a,
		    const Matrix<F, Dynamic, Dynamic>& A20,
		    const Matrix<F, Dynamic, Dynamic>& A11,
		    Matrix<F, Dynamic, Dynamic>* res) {
    *res = F(2) * ((a * a.transpose()).array() * A11.array());
    *res+=F(2) * ((A20 * a).array() * a.array()).matrix().asDiagonal();
  }

  template<class F>
  void Calc_ai_A_Bj_b(const Matrix<F, Dynamic, 1>& a1,
		      const Matrix<F, Dynamic, Dynamic>& A,
		      const Matrix<F, Dynamic, Dynamic>& B10,
		      const Matrix<F, Dynamic, 1>& b,
		      Matrix<F, Dynamic, Dynamic>* res) {
    *res =
      (a1 * (B10 * b).transpose()).array() * A.array() +
      (a1 * b.transpose()).array() * (A * B10.transpose()).array();
  }

  template<class F>
  void Calc_ai_b(const Matrix<F, Dynamic, 1>& a1,
		 const Matrix<F, Dynamic, 1>& b,
		 Matrix<F, Dynamic, 1>* res) {
    *res = a1.array() * b.array();
  }

  typedef std::complex<double> CD;
  template void Calc_a_Aj_b<double>(const VectorXd&,
				    const MatrixXd&,
				    const VectorXd&,
				    VectorXd*);  
  template void Calc_a_Aj_b<CD>(const VectorXcd&,
				const MatrixXcd&,
				const VectorXcd&,
				VectorXcd*);
  template void Calc_a_Ai_B_Aj_b<double>(const VectorXd&,
					 const MatrixXd&,
					 const MatrixXd&,
					 const VectorXd&,
					 MatrixXd*);
  template void Calc_a_Ai_B_Aj_b<CD>(const VectorXcd&,
				     const MatrixXcd&,
				     const MatrixXcd&,
				     const VectorXcd&,
				     MatrixXcd*);
  template void Calc_a_Aij_a<double>(const VectorXd&,
				     const MatrixXd&,
				     const MatrixXd&,
				     MatrixXd*);
  template void Calc_a_Aij_a<CD>(const VectorXcd&,
				 const MatrixXcd&,
				 const MatrixXcd&,
				 MatrixXcd*);
  template void Calc_ai_A_Bj_b<double>(const VectorXd&,
				       const MatrixXd&,
				       const MatrixXd&,
				       const VectorXd&,
				       MatrixXd*);
  template void Calc_ai_A_Bj_b<CD>(const VectorXcd&,
				   const MatrixXcd&,
				   const MatrixXcd&,
				   const VectorXcd&,
				   MatrixXcd*);
  template void Calc_ai_b<double>(const VectorXd&,
				  const VectorXd&,
				  VectorXd*);
  template void Calc_ai_b<CD>(const VectorXcd&,
			      const VectorXcd&,
			      VectorXcd*);  
}
