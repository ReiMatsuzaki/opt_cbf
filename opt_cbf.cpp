
#include <stdexcept>
#include <fstream>
#include <Eigen/LU>
#include <l2func.hpp>
#include <macros.hpp>
#include "l_algebra.hpp"
#include "opt_cbf.hpp"


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

  template<class BasisPrim, class DrivPrim>
  class OptCBF<BasisPrim, DrivPrim>::Impl {
  private:
    // ----------- type -------------------
    typedef typename BasisPrim::Field F;
    typedef typename vector<BasisPrim>::const_iterator IT;
    typedef LinearComb<BasisPrim> BasisLC;
    typedef LinearComb<DrivPrim>  DrivLC;
    typedef MatrixXcd M;
    typedef VectorXcd V;
  public:
    // ----------- Member Field ------------
    // inputs
    vector<BasisPrim> basis_set_;
    vector<BasisLC>   d_basis_set_;
    vector<BasisLC>   dd_basis_set_;
    Op<BasisPrim>     op_l_;
    DrivLC            driven_term_;
    // outputs
    VectorXcd coef_;
    
    // ----------- Constructors ------------
    Impl(const vector<BasisPrim>& us, const HLikeAtom<CD>& h_final,
	 const DrivLC& _driv, CD _energy) {
      basis_set_ = us;
      d_basis_set_.resize(us.size());
      dd_basis_set_.resize(us.size());
      op_l_ = h_final.HMinusEnergy<BasisPrim>(_energy);
      driven_term_ = _driv;

      VectorXcd zs;
      this->getZeta(&zs);
      this->updateZeta(zs);

    }
    // ---------- minor method -----------
    void getZeta(V* zs) {
      int num = basis_set_.size();
      zs->resize(num);
      for(int i = 0; i < num; i++) 
	(*zs)(i) = basis_set_[i].z();
    }
    void updateZeta(const V& zs) {

      int num = zs.rows();
      for(int i = 0; i < num; i++) {
	CD z = zs[i];
	BasisPrim ui(basis_set_[i].n(), z, Normalized);
	basis_set_[i]   = ui;
	d_basis_set_[i] = D1Normalized(ui);
	dd_basis_set_[i]= D2Normalized(ui);
      }
      
    }
    void computeMatrix(M* D00, M* D10, M* D20, M* D11) {

      int num = basis_set_.size();

      *D00 = MatrixXcd::Zero(num, num);
      *D10 = MatrixXcd::Zero(num, num);
      *D11 = MatrixXcd::Zero(num, num);
      *D20 = MatrixXcd::Zero(num, num);

      for(int j = 0; j < num; j++) {
	LinearComb<BasisPrim> l_uj   = op_l_(basis_set_[j]);
	LinearComb<BasisPrim> l_d_uj = op_l_(d_basis_set_[j]);

	(*D00)(j, j) = CIP(basis_set_[j], l_uj);
	(*D11)(j, j) = CIP(d_basis_set_[j], l_d_uj);

	for(int i = 0; i < num; i++) {
	  (*D10)(i, j) = CIP(d_basis_set_[i], l_uj);
	  (*D20)(i, j) = CIP(dd_basis_set_[i], l_uj);
	}

	for(int i = 0; i < j; i++) {
	  (*D00)(i, j) = (*D00)(j, i) = CIP(basis_set_[i], l_uj);
	  (*D11)(i, j) = (*D11)(j, i) = CIP(d_basis_set_[i], l_d_uj);
	}
      }
    }
    void computeVector(V* m0, V* m1, V* m2) {

      int num = basis_set_.size();
      
      *m0 = VectorXcd::Zero(num);
      *m1 = VectorXcd::Zero(num);
      *m2 = VectorXcd::Zero(num);

      for(int i = 0; i < num; i++) {
	(*m0)(i) =  CIP(basis_set_[i], driven_term_);
	(*m1)(i) =  CIP(d_basis_set_[i], driven_term_);
	(*m2)(i) =  CIP(dd_basis_set_[i], driven_term_);
      }

    }      
    // ---------- called from interface -------
    void Compute(const V& zs, CD* a, V* g, M* h) {

      // update orbital exponents (zs)
      this->updateZeta(zs);

      // prepare basic matrix and vector
      MatrixXcd D00, D10, D20, D11;
      VectorXcd m0, m1, m2;
      this->computeMatrix(&D00, &D10, &D20, &D11);
      this->computeVector(&m0, &m1, &m2);

      // compute coefficient
      coef_ = D00.fullPivLu().solve(m0);

      // compute gradient and hessian
      computeAlphaGradHess(coef_, 
			   D00, D10, D20, D11, m0, m1, m2, 
			   a, g, h);

    }
    BasisLC GetWaveFunction() const {

      LinearComb<BasisPrim> psi;
      int num_basis = basis_set_.size();
      for(int i = 0; i < num_basis; i++) 
	psi += coef_(i, 0) * basis_set_[i];
      return psi;      

    }
    void Display() const {
      cout << endl;
      cout << "========= OptCBF::Display ==========" << endl;
      cout << "Basis Set:" << basis_set_.size() << endl;
      for(int i = 0; i < basis_set_.size(); i++) 
	cout << "basis " << i << " : " << basis_set_[i] << endl;
      cout << "Driven Term:" << endl;
      for(int i= 0; i < driven_term_.size(); i++) {
	typename DrivPrim::Field c = driven_term_.coef_i(i);
	cout << "driv " << i <<  " : " << c;
	cout << " " << driven_term_.prim_i(i) << endl;
      }
      cout << "====================================" << endl;
    }
    void WritePsi(const string& fn, double rmax, double dr) const{

      ofstream ofs(fn.c_str());

      if(ofs.fail()) {
	string msg;
	SUB_LOCATION(msg);
	msg+="failed to open file: ";
	msg += string(fn);
	throw runtime_error(msg);
      }
      
      int num_grid(int(rmax/dr));
      LinearComb<BasisPrim> psi = this->GetWaveFunction();
      
      for(int i = 0; i < num_grid; i++) {

	double r = i * dr;
	CD     v = psi.at(r);
	ofs << r << ", " << v.real();
	ofs <<      ", " << v.imag() << endl;
	
      }

    }
    VectorXcd GetCoefs() const {
      return coef_;
    }
	 
  };

  template<class B, class D>
  OptCBF<B,D>::OptCBF(const vector<B>& a, const HLikeAtom<CD>& b,
		      const LinearComb<D>& c, CD d) :
    impl_(new Impl(a,b,c,d)) {}

  template<class B, class D> OptCBF<B,D>::~OptCBF() {
    delete impl_;
  }

  template<class B, class D>
  void OptCBF<B, D>::Compute(cV& zs, CD* a, VectorXcd* g, MatrixXcd* h) {

    impl_->Compute(zs, a, g, h);

    }
  template<class B, class D>
  LinearComb<B> OptCBF<B,D>::GetWaveFunction() const {
    return impl_->GetWaveFunction();
  }
  template<class B, class D> void OptCBF<B,D>::Display() const {
    impl_->Display();
  }
  template<class B, class D> void OptCBF<B,D>::WritePsi(string a, double b, double c) const {
    impl_->WritePsi(a, b, c);
  }
  template<class B, class D> VectorXcd OptCBF<B,D>::GetCoefs() const{
    return impl_->GetCoefs();
  }
  
  // ============= explicit instance ===============
  template class OptCBF<CSTO, CSTO>;
  template class OptCBF<CSTO, CGTO>;
  template class OptCBF<CGTO, CGTO>;
  template class OptCBF<CGTO, CSTO>;

  
}
