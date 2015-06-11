#include <iostream>
#include "driv.hpp"
#include <l2func.hpp>

using namespace std;
using namespace l2func;

namespace opt_cbf_h {

  template<class Basis>
  class HAtomPI<Basis>::Impl {
  private:
    // ------ type -----------
    typedef typename Basis::Field F;
    typedef typename LinearComb<Basis>::const_iterator IT;
    typedef LinearComb<Basis> LCBasis;

    // ------ field ----------
    int l_;   // angular quantum number;
    F   z_;   // charge;
    F   ene_; // energy;
    STOs driv_;
    
  public:
    // ------ constructors ---
    Impl(int _l, F _z, F _ene, STOs _driv ) :
      l_(_l), z_(_z), ene_(_ene), driv_(_driv) {
      IsPrimitive<Basis>();
      
    }
    ~Impl() { }

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
    F OpEle(const LCBasis& as, const Basis& b) {

      F acc(0);
      
      for(IT it = as.begin(), end_it = as.end(); it != end_it; ++it) {
	acc += it->first * this->OpEle(it->second, b);
      }
      return acc;
    }
    F OpEle(const LCBasis& as, const LCBasis& bs) {

      F acc(0);
      for(IT it = bs.begin(), it_end = bs.end(); it != it_end; ++it)
	acc += it->first * this->OpEle(as, it->second);
      return acc;
    }
    F DrivEle(const LCBasis& as) {

      F acc(0);
      for(IT it = as.begin(), end_it = as.end();
	  it != end_it; ++it) {
	acc += it->first * this->DrivEle(it->second);
      }
      return acc;
    }

    // ------- print --------------
    void Display() {
      cout << "l_z_e_d: " << l_ << ", " << z_ << ", " <<ene_ ;
      cout << ", " << driv_.size() << endl;
    }
  };

  template<class Basis>
  HAtomPI<Basis>::HAtomPI(int _l, F _z, F _ene, STOs _driv) :
    impl_(new typename HAtomPI<Basis>::Impl(_l, _z, _ene, _driv)) {}
  template<class Basis>
  HAtomPI<Basis>::~HAtomPI() { 
    delete impl_;
  }
  template<class Basis>
  typename Basis::Field HAtomPI<Basis>::OpEle
  (const Basis& a, const Basis& b) {
    return impl_->OpEle(a, b); 
  }
  template<class Basis>
  typename Basis::Field HAtomPI<Basis>::OpEle
  (const LinearComb<Basis>& a, const Basis& b) {
    return impl_->OpEle(a, b); 
  }
  template<class Basis>
  typename Basis::Field HAtomPI<Basis>::OpEle
  (const LinearComb<Basis>& a, const LinearComb<Basis>& b) {
    return impl_->OpEle(a, b); 
  }
  template<class Basis>
  typename Basis::Field HAtomPI<Basis>::DrivEle
  (const LinearComb<Basis>& a) {
    return impl_->DrivEle(a); 
  }
  template<class Basis>
  typename Basis::Field HAtomPI<Basis>::DrivEle
  (const Basis& a) {
    return impl_->DrivEle(a); 
  }
  template<class Basis>
  void HAtomPI<Basis>::Display() {
    impl_->Display();
  }
  

  // Explicit instance
  template class HAtomPI<RSTO>;
  template class HAtomPI<CSTO>;
  template class HAtomPI<RGTO>;
  template class HAtomPI<CGTO>;
    
}

