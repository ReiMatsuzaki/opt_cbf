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

namespace opt_cbf_h {

  // ------- forward declaration --------
  template<class Prim> class HAtomPI;

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
    virtual void Display() = 0;
    virtual void WritePsi(const string&,double,double)=0;
  };  

  // ====== calculator of alpha, gradient and hessian ==
  template<class Prim>
  class OptCBF : public IOptTarget {
  private:
    // ------ Member Field -------
    class Impl;
    Impl* impl_;

    // ------ uncopyable ---------
    OptCBF(const OptCBF&);
    
  public:
    // ------ constructor --------
    OptCBF(const vector<Prim>& _opt_basis_set,
	   HAtomPI<Prim>* _h_atom_pi);
    OptCBF(const vector<Prim>& _opt_basis_set,
	   const vector<Prim>& _fix_basis_set,
	   HAtomPI<Prim>* _h_atom_pi);
    ~OptCBF();

    // ------ method -----------
    // calculation of gradient and hessian from orbital exp.
    // this method will be used for optimization.
    void Compute(cV& zs, CD* a, VectorXcd* g, MatrixXcd* h);
    void computeMatrix (MatrixXcd* D00, MatrixXcd* D10,
			MatrixXcd* D20, MatrixXcd* D11);    
    void computeVector (VectorXcd* m0, VectorXcd* m1, 
			VectorXcd* m2);
    void Display();
    void WritePsi(const string&, double rmax, double dr);
    
  };
}

#endif
