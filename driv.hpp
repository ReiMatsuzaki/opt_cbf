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
  template<class Basis>
  class HAtomPI : public IDrivSystem<Basis> {
    
    // ------ type -----------
    typedef typename Basis::Field F;
    typedef LinearComb<ExpBasis<F, 1> > STOs;
    typedef typename LinearComb<Basis>::const_iterator IT;

    // ------ field ----------
    int l_;   // angular quantum number;
    F   z_;   // charge;
    F   ene_; // energy;
    STOs driv_;
    
  public:
    // ------ constructors ---
    HAtomPI(int _l, F _z, F _ene, STOs _driv ) :
      l_(_l), z_(_z), ene_(_ene), driv_(_driv) {}
    ~HAtomPI() { IsPrimitive<Basis>(); }

    // ------ matrix elements ------
    F OpEle(const Basis& a, const Basis& b) {

      F acc(0);
      acc += -0.5 * CIP(a, Op(OpDDr2<Basis>(), b));
      F ll = l_ * (l_ + 1) * 0.5;
      acc += ll    * CIP(a, Op(OpRm<Basis>(-2), b));
      acc += -z_   * CIP(a, Op(OpRm<Basis>(-1), b));
      acc += -ene_ * CIP(a, b);

      return acc;
    }
    F DrivEle(const Basis& a) {
      return CIP(a, driv_);
    }
    F OpEle(const LinearComb<Basis>& as, const Basis& b) {

      F acc(0);
      
      for(IT it = as.begin(), end_it = as.end(); it != end_it; ++it) {
	acc += it->first * this->OpEle(it->second, b);
      }
      return acc;
    }
    F OpEle(const LinearComb<Basis>& as, const LinearComb<Basis>& bs) {

      F acc(0);
      for(IT it = bs.begin(), it_end = bs.end(); it != it_end; ++it)
	acc += it->first * this->OpEle(as, it->second);
      return acc;
    }    F DrivEle(const LinearComb<Basis>& as) {

      F acc(0);
      for(IT it = as.begin(), end_it = as.end(); it != end_it; ++it)
	acc += it->first * this->DrivEle(it->second);
      return acc;
    }
  };
}

#endif
