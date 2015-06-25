#ifndef HATOM_HPP
#define HATOM_HPP

#include <stdexcept>
#include "prim.hpp"
#include "lcomb.hpp"

namespace l2func {

  template<class F=double>
  class HLikeAtom {
  private:
    // ---------- Field ------------
    int n_;
    F z_;
    int l_;

    // ---------- type -------------
    typedef LinearComb<ExpBasis<F,1> > STOs;
    
  public:
    // -------- constructor --------
    HLikeAtom(int _n, F _z, int _l);

    // -------- operator -----------
    template<class Prim>
    LinearComb<Prim> OperateHamiltonian(const Prim& o);

    template<class Prim>
    boost::function<LinearComb<Prim>(const Prim&)> Hamiltonian();

    // -------- state vector -------
    STOs EigenState();
    STOs DipoleInitLength(int l1);
    STOs DipoleInitVelocity(int l1);
    double EigenEnergy();
  };
}
#endif



