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
  using vector;
}

namespace opt_cbf_h {

  template<class Prim>
  class OptCBF {

    // ------ member field --------
  private:
    vector<Prim> basis_set_;
    HAtomPI<CD>  h_atom_pi_;
    
  public:
    // ------ constructor --------
    OptCBF(const vector<Prim>& _basis_set,
	   const HAtomPI<CD>&  _h_atom_pi) :
      basis_set_(_basis_set), h_atom_pi_(_h_atom_pi) {}

    // ------ method -----------
    void CalcGradHess(VectorXcd* g, MatrixXcd* h) {

      vector<Prim> d_bs;
      vector<Prim> dd_bs;

      calcDerivBasis(&d_bs, &dd_bs);
      
    }

    void calcDerivBasis(vector<Prim>* d_bs,
			vecetor<Prim>* dd_bs) {

    }
  }
}

#endif
