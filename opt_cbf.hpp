#ifndef OPT_CBF_HPP
#define OPT_CBF_HPP

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
}

namespace opt_cbf_h {

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
      it_fix_begin_(_opt_basis_set.begin()),
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
    void Compute(const vector<F>& zs, F* alpha,
		 VectorXcd* g, MatrixXcd* h) {

      // update orbital exponent (zeta)
      updateZeta(zs);

      // prepare matrix
      int num_all = basis_set_.size();
      int num_opt = numOptBasis();
      MatrixXcd D = MatrixXcd::Zero(num_all, num_all);
      VectorXcd m = VectorXcd::Zero(num_all);
      
      for(int i = 0; i < num_all; i++) {

	m(i, 0) = h_atom_pi_.DrivEle(basis_set_[i]);
	
	for(int j = 0; j < num_all; j++) {

	  D(i, j) = h_atom_pi_.OpEle(basis_set_[i], basis_set_[j]);
	  
	}
      }

      VectorXcd D_inv_m = D.fullPivLu().solve(m); 
      CD m_DInv_m = m.dot(D_inv_m);
      (*alpha) = m_DInv_m;
    }
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
    void updateZeta(const vector<F>& zs) {

      int num = zs.size();
      for(int i = 0; i < num; i++) {
	F z = zs[i];
	Prim f(basis_set_[i].n(), z, Normalized);
	basis_set_[i] = f;
	d_basis_set_[i]   = D1Normalized(f);
	dd_basis_set_[i]  = D2Normalized(f);
      }
    }
    void init() {
      
      // copy orbital exponent of opt_basis_set
      vector<F> zs;
      typedef typename vector<Prim>::const_iterator IT;
      for(IT it = basis_set_.begin(); it != it_fix_begin_; ++it) {
	zs.push_back(it->z());
	
      // call function to compute derivative basis set
      updateZeta(zs);
      
      }
    }
  };
}





#endif
