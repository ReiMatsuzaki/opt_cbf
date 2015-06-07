#ifndef DRIV_HPP
#define DRIV_HPP

// ======== forward decralation ============
namespace l2func {
  template<class F, int m> class ExpBasis;
  template<class Prim> class LinearComb;
}

namespace opt_cbf_h {

  // express driven type equation for hydrogen problem
  // (H-E)\psi = \phi
  template<class Basis>
  class IDrivSystem {
  private:
    typedef typename Basis::Field F;
    
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
    // ------- type -----------
    typedef typename Basis::Field F;
    typedef l2func::LinearComb<l2func::ExpBasis<F, 1> > STOs;
    typedef l2func::LinearComb<Basis> LCBasis;

    // ------- PIMPL idiom -----
    class Impl;
    Impl* impl_;

    // ------- uncopyable -------
    HAtomPI(const HAtomPI<Basis>&);
    HAtomPI& operator= (const HAtomPI<Basis>&);
    
  public:
    HAtomPI(int _l, F _z, F _ene, STOs _driv);
    
    ~HAtomPI();
    F OpEle(const Basis& a, const Basis& b);
    F OpEle(const LCBasis& a, const Basis& b);
    F OpEle(const LCBasis& a, const LCBasis& b);
    F DrivEle(const Basis& a);
    F DrivEle(const LCBasis& a);
    void Display();
  };
}

#endif
