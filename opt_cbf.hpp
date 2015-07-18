#ifndef OPT_CBF_HPP
#define OPT_CBF_HPP

#include <iostream>
#include <complex>
#include <vector>
#include <Eigen/Core>

// represents optimization target for CBF

namespace {
  using namespace Eigen;
  typedef std::complex<double> CD;
  using std::vector;
  using std::cout;
  using std::endl;
}

// ------- forward declaration --------
namespace l2func {
  template<class F> class HLikeAtom;
  template<class Prim> class LinearComb;
}


namespace opt_cbf_h {

  // ------- type -----------------------
  typedef const MatrixXcd cM;
  typedef const VectorXcd cV;

  // ------- logic ----------------------
  // compute mD^-1m and its gradient and hessian.
  void computeAlphaGradHess(cV& D_inv_m,
			    cM& D00, cM& D10, cM& D20, cM& D11,
			    cV& m0,  cV&m1,   cV& m2,
			    CD* a, VectorXcd* g, MatrixXcd* h);

  // ====== interface to access cSTO or cGTO ======
  class IOptTarget {
  public:
    virtual ~IOptTarget() {} ;
    virtual void Compute(const VectorXcd& zs, CD* a, VectorXcd* g, MatrixXcd* h) = 0;
    //    virtual void Display() = 0;
    //virtual void WritePsi(const string&,double,double)=0;
    //virtual VectorXcd GetCoefs() const = 0;
  };  

  // ====== calculator of alpha, gradient and hessian ==
  template<class BasisPrim, class DrivPrim>
  class OptCBF : public IOptTarget {
  private:
    // ------ Member Field -------
    class Impl;
    Impl* impl_;
    class ImplOld;

    // ------ uncopyable ---------
    OptCBF(const OptCBF&);
    
  public:
    // ------ constructor --------
    OptCBF(const vector<BasisPrim>& _opt_basis_set,
	   const l2func::HLikeAtom<CD>& h_final,
	   const l2func::LinearComb<DrivPrim>& _driv, 
	   CD _energy);
    ~OptCBF();

    // ------ method -----------
    // calculation of gradient and hessian from orbital exp.
    // this method will be used for optimization.
    void Compute(cV& zs, CD* a, VectorXcd* g, MatrixXcd* h);
    l2func::LinearComb<BasisPrim> GetWaveFunction() const;
    
  };
}

#endif
