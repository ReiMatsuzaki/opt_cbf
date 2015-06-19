
#include <stdexcept>
#include <fstream>
#include <Eigen/LU>
#include <l2func.hpp>
#include <macros.hpp>
#include "l_algebra.hpp"
#include "opt_cbf.hpp"
#include "driv.hpp"


using namespace l2func;
using std::ofstream;
using std::runtime_error;

namespace opt_cbf_h {
  
  void computeAlphaGradHess(cV& D_inv_m,
			    cM& D00, cM& D10, cM& D20, cM& D11,
			    cV& m0,  cV&m1,   cV& m2,
			    CD* a, 
			    VectorXcd* g, MatrixXcd* h) {
  
    //    VectorXcd D_inv_m = D00.fullPivLu().solve(m0);
    VectorXcd tmp = VectorXcd::Zero(m0.rows());

    // alpha
    *a = (m0.array() * D_inv_m.array()).sum();

    // grad
    Calc_a_Aj_b(D_inv_m, D10, D_inv_m, g);
    Calc_ai_b(m1, D_inv_m, &tmp);
    (*g) *= -1;
    (*g) += 2 * tmp;

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
    *h = tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6;
  }

  template<class Prim>
  class OptCBF<Prim>::Impl {
  private:
    // ------ type --------
    typedef typename Prim::Field F;
    typedef typename vector<Prim>::const_iterator IT;
    typedef LinearComb<Prim> LC;
    
  public:
    /* ------ member field --------
     * basis_set = {opt_basis1, ............... , opt_basisM }
     *           fix_basis1(= it_fix_begin), ...., fix_basisN,
    */
    vector<Prim>  basis_set_;
    IT            it_fix_begin_;
    HAtomPI<Prim>* h_atom_pi_;
    vector<LC>    d_basis_set_;    // derivative basis
    vector<LC>    dd_basis_set_;   // second derivative basis

    VectorXcd     coef_;

    // ------- constructor ---------
    Impl(const vector<Prim>& _opt_basis_set,
	 HAtomPI<Prim>* _h_atom_pi) :
      basis_set_(_opt_basis_set),
      it_fix_begin_(basis_set_.end()),
      h_atom_pi_(_h_atom_pi),
      d_basis_set_(_opt_basis_set.size()),
      dd_basis_set_(_opt_basis_set.size()) {

      this->init();
    }
    Impl(const vector<Prim>& _opt_basis_set,
	 const vector<Prim>& _fix_basis_set,
	 HAtomPI<Prim>*  _h_atom_pi) :
      h_atom_pi_(_h_atom_pi),
      d_basis_set_(_opt_basis_set.size()),
      dd_basis_set_(_opt_basis_set.size()) {
      
      // copy _opt_basis_set and
      basis_set_ = _opt_basis_set;
      it_fix_begin_ = basis_set_.end();
      basis_set_.insert(basis_set_.end(), 
			_fix_basis_set.begin(),
			_fix_basis_set.end());

      // initialize
      this->init();
      
    }
    ~Impl() {
      delete h_atom_pi_;
    }
    
    // ------- method ----------------
    void Compute(cV& zs, CD* a, VectorXcd* g, 
		 MatrixXcd* h) {

      // update orbital exponent (zeta)
      updateZeta(zs);

      // prepare matrix
      //      int num_all = basis_set_.size();
      //int num_opt = numOptBasis();

      MatrixXcd D00, D10, D20, D11;
      VectorXcd m0, m1, m2;
      //      VectorXcd tmp = VectorXcd::Zero(num_opt);
      

      // compute matrix and vector
      computeMatrix(&D00, &D10, &D20, &D11);
      computeVector(&m0, &m1, &m2);
      coef_ = D00.fullPivLu().solve(m0);

      computeAlphaGradHess(coef_,
			   D00, D10, D20, D11, m0, m1, m2, 
			   a, g, h);
      
    }
    void computeMatrix
    (MatrixXcd* D00, MatrixXcd* D10,
     MatrixXcd* D20, MatrixXcd* D11) {

      // typedef typename vector<Prim>::const_iterator IT;

      int num_all = basis_set_.size();
      int num_opt = numOptBasis();

      *D00 = MatrixXcd::Zero(num_all, num_all);
      *D10 = MatrixXcd::Zero(num_opt, num_all);
      *D11 = MatrixXcd::Zero(num_opt, num_opt);
      *D20 = MatrixXcd::Zero(num_opt, num_all);

      for(int i = 0; i < num_all; i++) {
	CD d_00_ii = h_atom_pi_->OpEle(basis_set_[i], 
				    basis_set_[i]);
	CD d_11_ii = h_atom_pi_->OpEle(d_basis_set_[i],
				       d_basis_set_[i]);
	(*D00)(i, i) = d_00_ii;
	(*D11)(i, i) = d_11_ii;
	for(int j = i+1; j < num_all; j++) {
	  CD d_00_ij = h_atom_pi_->OpEle(basis_set_[i], 
					 basis_set_[j]);
	  CD d_11_ij = h_atom_pi_->OpEle(d_basis_set_[i], 
					 d_basis_set_[j]);
	  (*D00)(i, j) = d_00_ij;
	  (*D00)(j, i) = d_00_ij;
	  (*D11)(i, j) = d_11_ij;
	  (*D11)(j, i) = d_11_ij;
	}
      }

      for(int i = 0; i < num_opt; i++)
	for(int j = 0; j < num_all; j++) {
	  (*D10)(i,j) = h_atom_pi_->OpEle(d_basis_set_[i],  
					  basis_set_[j]);
	  (*D20)(i,j) = h_atom_pi_->OpEle(dd_basis_set_[i], 
					  basis_set_[j]);
	}
      
    }  
    void computeVector
    (VectorXcd* m0, VectorXcd* m1, VectorXcd* m2) {

      int num_all = basis_set_.size();
      int num_opt = numOptBasis();
      
      *m0 = VectorXcd::Zero(num_all);
      *m1 = VectorXcd::Zero(num_opt);
      *m2 = VectorXcd::Zero(num_opt);

      for(int i = 0; i < num_all; i++)
	(*m0)(i, 0) = h_atom_pi_->DrivEle(basis_set_[i]);
      for(int i = 0; i < num_opt; i++) {
	(*m1)(i, 0) = h_atom_pi_->DrivEle(d_basis_set_[i]);
	(*m2)(i, 0) = h_atom_pi_->DrivEle(dd_basis_set_[i]);
      }
      
    }
    int numOptBasis()  {
    
      typedef typename vector<Prim>::const_iterator IT;
      return std::distance(static_cast<IT>(basis_set_.begin()),
			   it_fix_begin_);
    }
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
    void init() {
      
      // copy orbital exponent of opt_basis_set
      VectorXcd zs = VectorXcd::Zero(basis_set_.size());

      int i = 0;
      typedef typename vector<Prim>::const_iterator IT;
      for(IT it = basis_set_.begin(); it != it_fix_begin_; 
	  ++it, ++i) 
	zs(i) = it->z();
      
      // call function to compute derivative basis set
      updateZeta(zs);
    }
    void Display() {
      cout << "OptCBF_Display" << endl;
      h_atom_pi_->Display();
    }
    LinearComb<Prim> WaveFunc() const {
      
      LinearComb<Prim> psi;
      int num_basis = basis_set_.size();
      for(int i = 0; i < num_basis; i++) 
	psi += coef_(i, 0) * basis_set_[i];
      return psi;
      
    }
    void WritePsi(const string& fn, double rmax, double dr) {

      ofstream ofs(fn.c_str());

      if(ofs.fail()) {
	string msg;
	SUB_LOCATION(msg);
	msg+="failed to open file: ";
	msg += string(fn);
	throw runtime_error(msg);
      }
      
      int num_grid(int(rmax/dr));
      LinearComb<Prim> psi = this->WaveFunc();
      
      for(int i = 0; i < num_grid; i++) {

	double r = i * dr;
	CD     v = AtX(r, psi);
	ofs << r << ", " << v.real();
	ofs <<      ", " << v.imag() << endl;
	
      }
    }
    VectorXcd GetCoefs() const {
      return coef_;
    }
  };


  template<class Prim>
  OptCBF<Prim>::OptCBF(const vector<Prim>& _opt_basis_set,
		       HAtomPI<Prim>* _h_atom_pi) :
    impl_(new Impl(_opt_basis_set, _h_atom_pi)) {}

  template<class Prim>
  OptCBF<Prim>::OptCBF(const vector<Prim>& _o,
		       const vector<Prim>& _f,
		       HAtomPI<Prim>* _h) :
    impl_(new Impl(_o, _f, _h)) {}

  template<class Prim>
  OptCBF<Prim>::~OptCBF() {
    delete impl_;
  }

  template<class Prim>
  void OptCBF<Prim>::Compute(cV& zs, CD* a, VectorXcd* g, 
			     MatrixXcd* h) {

    impl_->Compute(zs, a, g, h);

    }
  template<class Prim>
  void OptCBF<Prim>::computeMatrix
  (MatrixXcd* D00, MatrixXcd* D10,
   MatrixXcd* D20, MatrixXcd* D11) {
    impl_->computeMatrix(D00, D10, D20, D11);
  }
  template<class Prim>
  void OptCBF<Prim>::computeVector 
    (VectorXcd* m0, VectorXcd* m1, VectorXcd* m2) {
      impl_->computeVector(m0, m1, m2);
    }

  template<class Prim>
  void OptCBF<Prim>::Display() {
    impl_->Display();
  }
  
  template<class Prim>
  void OptCBF<Prim>::WritePsi(const string& fn, double rm, double dr) {
    impl_->WritePsi(fn, rm, dr);
  }
  
  template<class Prim>
  VectorXcd OptCBF<Prim>::GetCoefs() const {
    return impl_->GetCoefs();
  }
  
  // ============= explicit instance ===============
  template class OptCBF<CSTO>;
  template class OptCBF<CGTO>;
  
}
