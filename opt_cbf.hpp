#ifndef OPT_CBF_HPP
#define OPT_CBF_HPP

#include <iostream>
#include <complex>
#include <vector>
#include <Eigen/Core>

// WARNING!! This is old version of opt_target.hpp
// Look opt_target instead of this file.
// represents optimization target for CBF


namespace {
  using namespace Eigen;
  typedef std::complex<double> CD;
  using std::vector;
  using std::cout;
  using std::endl;
  using std::string;
}

// ------- forward declaration --------
namespace l2func {
  template<class F> class HLikeAtom;
  template<class Prim> class LinearComb;
}


namespace opt_cbf_h {

  // ====== interface to access cSTO or cGTO ======
  class IOptTarget {
  public:
    virtual ~IOptTarget() {} ;
    virtual void Compute(const VectorXcd& zs, CD* a, VectorXcd* g, MatrixXcd* h) = 0;
    virtual void Display() const = 0;
    virtual void WritePsi(string,double,double) const = 0;
    virtual VectorXcd GetCoefs() const = 0;
  };  

  // ====== calculator of alpha, gradient and hessian ==
  template<class BasisPrim, class DrivPrim> class OptCBF : public IOptTarget {
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
    void Compute(VectorXcd& zs, CD* a, VectorXcd* g, MatrixXcd* h);
    l2func::LinearComb<BasisPrim> GetWaveFunction() const;
    void Display() const;
    void WritePsi(string, double, double) const;
    VectorXcd GetCoefs() const;
  };
}

#endif
