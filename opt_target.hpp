#ifndef OPT_TARGET_TEMPLATE_H
#define OPT_TARGET_TEMPLATE_H

#include <vector>
#include <Eigen/Core>
#include <l2func.hpp>

namespace {
  typedef Eigen::MatrixXcd M;
  typedef Eigen::VectorXcd V;  
  typedef const Eigen::MatrixXcd cM;
  typedef const Eigen::VectorXcd cV;  
  typedef std::complex<double> CD;  
  using std::string;
}

namespace opt_cbf_h {

  // ==== external function ====
  // compute vD^-1u and its gradient and hessian.
  void AlphaGradHess(cV& D_inv_m,
		     cM& D00, cM& D10, cM& D01, 
		     cM& D02, cM& D20, cM& D11,
		     cV& v0,  cV& v1,  cV& v2,
		     cV& u0,  cV& u1,  cV& u2,
		     CD* a, V* g, M* h);
  
  // ==== Interface to access any basis ====
  class IOptTarget {
  public:
    virtual ~IOptTarget() {}
    virtual void Compute(cV& zs, CD* a, V* g, M* h) = 0;
    //    virtual void Display() const = 0;
    //    virtual void WritePsi(string,double,double) = 0;
    //    virtual V    GetCoefs() const = 0;
  };
  
  // ==== Specific Basis ====
  /*
    FuncR : Primitive function to represent R
    FuncB : Primitive function for basis 
    FuncS : Primitive function for represent source term S
    OpL   : oprator for linear operaotr L
    
   */
  template<class FuncR, class FuncB, class FuncS, class OpL> 
  class OptTarget :public IOptTarget{
  private:
    // ---- type ----
    //    typedef typename FuncB::Field F;
    typedef typename std::vector<FuncB> FuncBs;
    //    typedef typename FuncBs::const_iterator IT;

    // ---- Implement Pointer ----
    class Impl;
    Impl *impl_;
    
  public:
    // ---- Constructors ----
    OptTarget(const FuncR& _r_func, 
	      const FuncBs& _basis_set,
	      const FuncS& _s_func,
	      const OpL& _op_l);
    ~OptTarget();

    // ---- Method ----
    virtual void Compute(cV& zs, CD* a, V* g, M* h);
    void Display();
    void WritePsi(string,double,double);
    V    GetCoefs() const;
  };

  // ==== type ====  
  typedef OptTarget<l2func::HLengthType<1,0,1,CD>::type,
		    l2func::NormalCSTO,
		    l2func::HLengthType<1,0,1,CD>::type,
		    l2func::HminusEOp<1,CD>::type> OptSTOLength;

}

#endif
