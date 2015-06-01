#ifndef OPT_CBF_HPP
#define OPT_CBF_HPP

#include <iostream>
#include <complex>
#include <vector>
#include <Eigen/Core>
#include <l2func.hpp>
#include "driv.hpp"

// represents optimization target for CBF

namespace {
  using namespace Eigen;
  typedef std::complex<double> CD;
  using std::vector;
  using std::cout;
  using std::endl;
}

namespace opt_cbf_h {

  typedef const MatrixXcd cM;
  typedef const VectorXcd cV;

  // (it is not used)
  // difference between matrix and its transpose
  void writeDiffMatrixAndItsTranspose(cM& m) {

    cout << (m - m.transpose()).array().abs().sum() << endl;
    
  }

  // compute mD^-1m and its gradient and hessian.
  void computeAlphaGradHess(cM& D00, cM& D10, cM& D20, cM& D11,
			    cV& m0,  cV&m1,   cV& m2,
			    CD* alpha, VectorXcd* grad, MatrixXcd* hess) {
  
    VectorXcd D_inv_m = D00.fullPivLu().solve(m0);
    VectorXcd tmp = VectorXcd::Zero(m0.rows());

    // alpha
    *alpha = (m0.array() * D_inv_m.array()).sum();

    // grad
    Calc_a_Aj_b(D_inv_m, D10, D_inv_m, grad);
    Calc_ai_b(m1, D_inv_m, &tmp);
    (*grad) *= -1;
    (*grad) += 2 * tmp;

    // hess
    MatrixXcd Dinv = D00.inverse();
    MatrixXcd tmp1 = (2*m2.array()*D_inv_m.array()).matrix().asDiagonal();

    MatrixXcd tmp2 = 2 * (m1 * m1.transpose()).array() * Dinv.array();
    MatrixXcd tmp3;
    Calc_a_Aij_a(D_inv_m, D20, D11, &tmp3);
    tmp3 *= -1;
    MatrixXcd tmp4;
    Calc_ai_A_Bj_b(m1, Dinv, D10, D_inv_m, &tmp4);
    tmp4 *= -2;
    MatrixXcd tmp5;
    tmp5 = tmp4.transpose();
    MatrixXcd tmp6;
    Calc_a_Ai_B_Aj_b(D_inv_m, D10, Dinv, D_inv_m, &tmp6);
    tmp6 *= 2;
    *hess = tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6;
  }

  template<class Prim>
  class OptCBF {

    // ------ type --------
    typedef typename Prim::Field F;
    typedef LinearComb<Prim> LC;

    // ------ member field --------
    // basis_set = {opt_basis1, ............... , opt_basisM }
    //              fix_basis1(= it_fix_begin), ...., fix_basisN,
  private:
    vector<Prim> basis_set_; // basis set which will be optimized.
    typename vector<Prim>::const_iterator it_fix_begin_;
    HAtomPI<Prim> h_atom_pi_;
    vector<LC> d_basis_set_;    // derivative basis
    vector<LC> dd_basis_set_;   // second derivative basis
    
  public:
    // ------ constructor --------
    OptCBF(const vector<Prim>& _opt_basis_set,
	   const HAtomPI<Prim>& _h_atom_pi) :
      basis_set_(_opt_basis_set),
      it_fix_begin_(basis_set_.end()),
      h_atom_pi_(_h_atom_pi),
      d_basis_set_(_opt_basis_set.size()),
      dd_basis_set_(_opt_basis_set.size()) {

      init();
      
    }
    OptCBF(const vector<Prim>& _opt_basis_set,
	   const vector<Prim>& _fix_basis_set,
	   const HAtomPI<CD>&  _h_atom_pi) :
    h_atom_pi_(_h_atom_pi),
    d_basis_set_(_opt_basis_set.size()),
    dd_basis_set_(_opt_basis_set.size()) {
      
      // copy _opt_basis_set and
      basis_set_ = _opt_basis_set;
      it_fix_begin_ = basis_set_.end();
      basis_set_.insert(basis_set_.end(), _fix_basis_set.begin(),
			_fix_basis_set.end());

      // initialize
      init();
      
    }

    // ------ method -----------
    // calculation of gradient and hessian from orbital exponent.
    // this method will be used for optimization.
    void Compute(const VectorXcd& zs, F* a, VectorXcd* g, MatrixXcd* h) {

      // update orbital exponent (zeta)
      updateZeta(zs);

      // prepare matrix
      //      int num_all = basis_set_.size();
      int num_opt = numOptBasis();

      MatrixXcd D00, D10, D20, D11;
      VectorXcd m0, m1, m2;
      VectorXcd tmp = VectorXcd::Zero(num_opt);

      // compute matrix and vector
      computeMatrix(&D00, &D10, &D20, &D11);
      computeVector(&m0, &m1, &m2);

      computeAlphaGradHess(D00, D10, D20, D11, m0, m1, m2, a, g, h);
      
    }
    // compute matrix and vector
    void computeMatrix(MatrixXcd* D00, MatrixXcd* D10,
		       MatrixXcd* D20, MatrixXcd* D11) {

      int num_all = basis_set_.size();
      int num_opt = numOptBasis();

      *D00 = MatrixXcd::Zero(num_all, num_all);
      *D10 = MatrixXcd::Zero(num_opt, num_all);
      *D11 = MatrixXcd::Zero(num_opt, num_opt);
      *D20 = MatrixXcd::Zero(num_opt, num_all);

      for(int i = 0; i < num_all; i++) 
	for(int j = 0; j < num_all; j++) 
	  (*D00)(i, j) = h_atom_pi_.OpEle(basis_set_[i], basis_set_[j]);

      for(int i = 0; i < num_opt; i++)
	for(int j = 0; j < num_opt; j++)
	  (*D11)(i,j) = h_atom_pi_.OpEle(d_basis_set_[i],d_basis_set_[j]);
      for(int i = 0; i < num_opt; i++)
	for(int j = 0; j < num_all; j++) {
	  (*D10)(i,j) = h_atom_pi_.OpEle(d_basis_set_[i],  basis_set_[j]);
	  (*D20)(i,j) = h_atom_pi_.OpEle(dd_basis_set_[i], basis_set_[j]);
	}
    }
    void computeVector(VectorXcd* m0, VectorXcd* m1, VectorXcd* m2) {

      int num_all = basis_set_.size();
      int num_opt = numOptBasis();
      
      *m0 = VectorXcd::Zero(num_all);
      *m1 = VectorXcd::Zero(num_opt);
      *m2 = VectorXcd::Zero(num_opt);

      for(int i = 0; i < num_all; i++)
	(*m0)(i, 0) = h_atom_pi_.DrivEle(basis_set_[i]);
      for(int i = 0; i < num_opt; i++) {
	(*m1)(i, 0) = h_atom_pi_.DrivEle(d_basis_set_[i]);
	(*m2)(i, 0) = h_atom_pi_.DrivEle(dd_basis_set_[i]);
      }
      
    }
    // number of basis for optimization
    int numOptBasis() {
      /*
      int acc(0);
      for(vector<Prim>::const_iterator it = basis_set_.begin();
	  it != it_fix_begin_; ++it)
	acc++;
      return acc;
      */
      typedef typename vector<Prim>::const_iterator IT;
      return std::distance(static_cast<IT>(basis_set_.begin()),
			   it_fix_begin_);
    }
    // update orbital exponents
    void updateZeta(const VectorXcd& zs) {

      int num = zs.rows();
      for(int i = 0; i < num; i++) {
	CD z = zs[i];
	Prim ui(basis_set_[i].n(), z, Normalized);
	basis_set_[i]   = ui;
	d_basis_set_[i] = D1Normalized(ui);
	dd_basis_set_[i]= D2Normalized(ui);
      }
      
    }
    /*
    void updateZeta(const vector<F>& zs) {

      VectorXcd vec_zs = Map<VectorXcd>(&zs[0], zs.size());

      // redirection
      this->updateZeta(vec_zs);
      
      int num = zs.size();
      for(int i = 0; i < num; i++) {
	F z = zs[i];
	Prim f(basis_set_[i].n(), z, Normalized);
	basis_set_[i] = f;
	d_basis_set_[i]   = D1Normalized(f);
	dd_basis_set_[i]  = D2Normalized(f);
      }

    }
    */
    // initialize
    void init() {
      
      // copy orbital exponent of opt_basis_set
      VectorXcd zs = VectorXcd::Zero(basis_set_.size());
      //      vector<F> zs;
      int i = 0;
      typedef typename vector<Prim>::const_iterator IT;
      for(IT it = basis_set_.begin(); it != it_fix_begin_; ++it, ++i) 
	zs(i) = it->z();
      
      // call function to compute derivative basis set
      updateZeta(zs);
    }
  };
}

#endif
