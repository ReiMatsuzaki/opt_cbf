#ifndef DRIV_HPP
#define DRIV_HPP

#include <vector>
#include <l2func.hpp>

namespace {
  using std::vector;
  using namespace l2func;
}

namespace opt_cbf_h {

  // express driven type equation for hydrogen problem
  // (H-E)\psi = \phi
  
  template<class F>
  class HAtomPI {
    typedef LinearComb<ExpBasis<F, 1> > STOs;
    
    int l_;   // angular quantum number;
    F   z_;   // charge;
    F   ene_; // energy;
    LinearComb<ExpBasis<F, 1> > driv_;
    
  public:
    HAtomPI(int _l, F _z, F _ene, STOs _driv ) :
      l_(_l), z_(_z), ene_(_ene), driv_(_driv) {}
    template<class Basis>
    F OpEle(Basis& a, Basis& b) {

      F acc(0);
      acc += -0.5 * CIP(a, Op(OpDDr2<RSTO>(), b));
      F ll = l_ * (l_ + 1) * 0.5;
      acc += ll    * CIP(a, Op(OpRm<RSTO>(-2), b));
      acc += -z_   * CIP(a, Op(OpRm<RSTO>(-1), b));
      acc += -ene_ * CIP(a, b);

      return acc;
    }
    template<class Basis>
    F DrivEle(Basis& a) {
      return CIP(a, driv_);
    }
  };
}

#endif
