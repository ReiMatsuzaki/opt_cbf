#include <Eigen/LU>
#include "opt_target.hpp"
#include "l_algebra.hpp"

namespace opt_cbf_h {

  void AlphaGrad(cV& d, cV& c,
		 cM& D10,
		 cV& r1,
		 cV& s0, cV& s1,
		 CD *a, V *g) {

    // alpha
    *a = (d.array() * s0.array()).sum();

    // alpha
    *a = (d.array() * s0.array()).sum();

    // grad
    V g1 = Calc_a_Aj_b(d, D10, c);
    V g2 = Calc_ai_b(r1, c);
    V g3 = Calc_ai_b(s1, d);
      
    *g = g2 + g3 - g1;
  }

  void AlphaGradHess(
		     cM& D00, cM& D10,
		     cM& D20, cM& D11,
		     cV& r0,  cV& r1,  cV& r2,
		     cV& s0,  cV& s1,  cV& s2,
		     CD* a, V* g, M* h) {

    // basic
    M U = D00.inverse();
    V c = U*s0;
    V d = U*r0;

    // alpha
    *a = (d.array() * s0.array()).sum();

    // grad
    V g1 = Calc_a_Aj_b(d, D10, c);
    V g2 = Calc_ai_b(r1, c);
    V g3 = Calc_ai_b(s1, d);
      
    *g = g2 + g3 - g1;

    // hess
    M rij_c       = (r2.array() * c.array()).matrix().asDiagonal();
    M ri_U_Lj_c   = Calc_ai_A_Bj_b(r1, U, D10, c);
    M ri_U_sj     = (r1 * s1.transpose()).array() * U.array();
    M d_Li_U_Lj_c = Calc_a_Ai_B_Aj_b(d, D10, U, c);
    M d_Lij_c     = Calc_a_Aij_b(d, D20, D11, c);
    M d_Li_U_sj   = Calc_ai_A_Bj_b(s1, U, D10, d);
    M d_sij       = (s2.array() * d.array()).matrix().asDiagonal();

    *h = rij_c - d_Lij_c + d_sij + 2*(-ri_U_Lj_c + ri_U_sj -ri_U_Lj_c + d_Li_U_Lj_c );

  }

  template<class FuncR, class FuncB, class FuncS, class OpL> 
  class OptTarget<FuncR, FuncB, FuncS, OpL>::Impl {
  private:
    
  public:
    // ---- Input Field ----
    FuncR   r_func_;
    FuncBs  basis_set_;
    FuncS   s_func_;
    OpL     op_l_;
    
    typedef typename FuncB::FuncDerivOne DFuncB;
    typedef typename FuncB::FuncDerivTwo DDFuncB;
    std::vector<DFuncB >  d_basis_set_;
    std::vector<DDFuncB>  dd_basis_set_;

    // ---- Output Field ----
    V coef_, coef_star_;

    // ---- Constructors ----
    Impl(const FuncR& _r_func, 
	 const FuncBs& _basis_set,
	 const FuncS& _s_func,
	 const OpL& _op_l
	 ) : r_func_(_r_func), 
	     basis_set_(_basis_set), 
	     s_func_(_s_func),
	     op_l_(_op_l) {

      V zs;
      this->getZeta(&zs);
      this->updateZeta(zs);

    }
    void initDBasis() {
      
      for(typename FuncBs::const_iterator it = basis_set_.begin();
	  it != basis_set_.end(); ++it) {

	d_basis_set_.push_back(it->DerivParamOne());
	dd_basis_set_.push_back(it->DerivParamTwo());

      }
    }
    void getZeta(V* zs) {

      int num = this->basis_set_.size();
      zs->resize(num);
      for(int i = 0; i < num; i++) {
	(*zs)(i) = basis_set_[i].z();
      }

    }
    void updateZeta(const V& zs) {

      int num = zs.rows();
      for(int i = 0; i < num; i++) {
	FuncB ui(this->basis_set_[i]);
	ui.set_z(zs[i]);
	
	this->basis_set_[i] = ui;
	this->d_basis_set_[i] = ui.DerivParamOne();
	this->dd_basis_set_[i] = ui.DerivParamTwo();
      }
    }
    void computeMatrix(M* D00, M* D10, M* D20, M* D11) {

      int num = basis_set_.size();

      *D00 = Eigen::MatrixXcd::Zero(num, num);
      *D10 = Eigen::MatrixXcd::Zero(num, num);
      *D20 = Eigen::MatrixXcd::Zero(num, num);
      *D11 = Eigen::MatrixXcd::Zero(num, num);

      for(int i = 0; i < num; i++) {
	for(int j = 0; j < num; j++) {

	  (*D00)(i ,j) = CIP(basis_set_[i], op_l_, basis_set_[j]);
	  (*D10)(i ,j) = CIP(d_basis_set_[i], op_l_, basis_set_[j]);
	  (*D20)(i ,j) = CIP(dd_basis_set_[i], op_l_, basis_set_[j]);
	  (*D11)(i ,j) = CIP(d_basis_set_[i], op_l_, d_basis_set_[j]);

	}
      }
      

    }
    void computeVector(V* r0, V* r1, V* r2, V* s0, V* s1, V* s2) {
      
      int num = basis_set_.size();
            
      *r0 = Eigen::VectorXcd::Zero(num);
      *r1 = Eigen::VectorXcd::Zero(num);
      *r2 = Eigen::VectorXcd::Zero(num);
      *s0 = Eigen::VectorXcd::Zero(num);
      *s1 = Eigen::VectorXcd::Zero(num);
      *s2 = Eigen::VectorXcd::Zero(num);

      for(int i = 0; i < (int)basis_set_.size(); i++) {

	(*r0)(i) = CIP(r_func_, basis_set_[i]);
	(*r1)(i) = CIP(r_func_, d_basis_set_[i]);
	(*r2)(i) = CIP(r_func_, dd_basis_set_[i]);
	(*s0)(i) = CIP(r_func_, basis_set_[i]);
	(*s1)(i) = CIP(d_basis_set_[i], s_func_);
	(*s2)(i) = CIP(dd_basis_set_[i], s_func_);

      }



    }
    void Compute(const V& zs, CD* a, V* g, M* h) {

      // update
      this->updateZeta(zs);
      M D00, D10, D20, D11;
      V r0, r1, r2, s0, s1, s2;

      // compute part
      this->computeMatrix(&D00, &D10, &D20, &D11);
      this->computeVector(&r0, &r1, &r2, &s0, &s1, &s2);

      CD aa;
      V gg;
      M hh;
      AlphaGradHess(
		    D00, D10, D20, D11,
		    r0, r1, r2, s0, s1, s2,
		    &aa, &gg, &hh);
      *a = aa;
      *g = gg;
      *h = hh;
    }
  };

  // ---- Method ----
  template<class FuncR, class FuncB, class FuncS, class OpL> 
  OptTarget<FuncR, FuncB, FuncS, OpL>::OptTarget(const FuncR& _r_func, 
						 const FuncBs& _basis_set,
						 const FuncS& _s_func,
						 const OpL& _op_l):
    impl_(new Impl(_r_func, _basis_set, _s_func, _op_l)) {}
  template<class FuncR, class FuncB, class FuncS, class OpL> 
  OptTarget<FuncR, FuncB, FuncS, OpL>::~OptTarget() { delete impl_; }

  template<class FuncR, class FuncB, class FuncS, class OpL> 
  void OptTarget<FuncR, FuncB, FuncS, OpL>::Compute(cV& zs, CD* a, V* g, M* h) {
    impl_->Compute(zs, a, g, h);
  }
  template<class FuncR, class FuncB, class FuncS, class OpL> 
  void OptTarget<FuncR, FuncB, FuncS, OpL>::Display() {
    
  }
  template<class FuncR, class FuncB, class FuncS, class OpL> 
  void OptTarget<FuncR, FuncB, FuncS, OpL>::WritePsi(string,double,double) {

  }
  template<class FuncR, class FuncB, class FuncS, class OpL> 
  V OptTarget<FuncR, FuncB, FuncS, OpL>::GetCoefs() const {
    return V::Zero(2);
  }

  // ==== Explicit Decraltion ====
  typedef std::complex<double> CD;
  template
  class OptTarget<l2func::HLengthType<1,0,1,CD>::type,
		  l2func::NormalCSTO,
		  l2func::HLengthType<1,0,1,CD>::type,
		  l2func::HminusEOp<0,CD>::type>;
		  
}
