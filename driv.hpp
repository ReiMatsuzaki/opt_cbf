#ifndef DRIV_HPP
#define DRIV_HPP

#include <vector>
#include <l2func.hpp>

//namespace {
using std::vector;
using namespace l2func;
//}

namespace opt_cbf_h {

  // express driven type equation for hydrogen problem
  // (H-E)\psi = \phi
  template<class Basis>
  class IDrivSystem {
  private:
    typedef typename Basis::Field F;
    typedef typename LinearComb<Basis>::const_iterator IT;
    
  public:
    virtual ~IDrivSystem() {}
    virtual F OpEle(const Basis& a, const Basis& b) = 0;
    virtual F DrivEle(const Basis& a) = 0;
  };
  
  // driven type equation for radial part of 
  // Hydrogen atom photoionization.
  // PIMPL idiom is used.
  template<class Basis>
  class HAtomPI : public IDrivSystem<Basis> {
  private:
    typedef LinearComb<ExpBasis<F, 1> > STOs;
    class Impl;
    Impl* impl_;
    
  public:
    HAtomPI(int _l, F _z, F _ene, STOs _driv);
    ~HAtomPI();
  };

}

#endif
