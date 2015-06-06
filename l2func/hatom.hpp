#ifndef HATOM_HPP
#define HATOM_HPP

#include <stdexcept>
#include "prim.hpp"
#include "lcomb.hpp"

namespace l2func {

  template<class F=double>
  class HydrogenLikeAtom {
  private:
    F z_;
    int l_;
  public:
    HydrogenLikeAtom(F _z, int _l) :
      z_(_z), l_(_l) {}
  };

  // ----------- Hamiltonian -------------------
  template<class Prim>
  LinearComb<Prim>
  OperateHAtomRadial(typename Prim::Field z, int l, const Prim& o) {
      IsPrimitive<Prim>();
      LinearComb<Prim> res;
      res += Op(OpCst<Prim>(-0.5), Op(OpDDr2<Prim>(), o));
      res += Op(OpCst<Prim>(-z), Op(OpRm<Prim>(-1), o));
      if(l != 0)
	res += Op(OpCst<Prim>(l*(l+1)/2.0), Op(OpRm<Prim>(-2), o));
      return res;        
  }
  template<class Prim>
  function<LinearComb<Prim>(const Prim&)>
  OpHAtomRadial(typename Prim::Field z, int l) {
    return bind(OperateHAtomRadial<Prim>, z, l, _1);
  }  

  // ----------- eigen state and dipole init ----
  template<class F>
  LinearComb<ExpBasis<F,1> > HAtomEigenState(int n, int l) {
    LinearComb<ExpBasis<F, 1> > lc;

    if( n == 1 && l==0) {
      
      lc += ExpBasis<F,1>(2.0, 1, 1.0);
      
    } else if ( n== 2 && l == 0) {
      
      lc += ExpBasis<F,1>(1.0/sqrt(2.0), 1, 0.5);
      lc += ExpBasis<F,1>(-1.0/(2.0*sqrt(2.0)), 2, 0.5);
      
    } else if (n == 2 && l == 1) {
      
      lc += ExpBasis<F,1>(1.0/(2.0 * sqrt(6.0)), 2, 0.5);

    } else if (n == 3 && l == 0) {

      lc += ExpBasis<F,1>(4.0 / (81.0 * sqrt(30.0)), 3, 1.0/3.0);
      
    } else if (n == 3 && l == 1) {
      double c = 8.0 / (27.0 * sqrt(6.0));
      double z = 1.0 / 3.0;
      lc += ExpBasis<F,1>(c,      2, z);
      lc += ExpBasis<F,1>(-c/6.0, 3, z);
      
    } else {

      string msg;
      msg = "inputted n and l is not supported for HAtomEigenState";
      throw std::invalid_argument(msg);

    }

    return lc;
  }
  template<class F>
  LinearComb<ExpBasis<F,1> > HAtomDipoleInitL(int n0,int l0,int l1) {

    if(l0 != l1 + 1 && l0 != l1 - 1) {

      string msg;
      msg = "|l0 - l1| != 1 in HAtomDipoleInitL";
      throw std::invalid_argument(msg);

    }

    typedef ExpBasis<F, 1> Prim;
    typedef LinearComb<Prim> LC;

    LC psi_n = HAtomEigenState<F>(n0, l0);
    LC res = Op(OpRm<Prim>(1), psi_n);
    return res;
    
  }
  template<class F>
  LinearComb<ExpBasis<F,1> > HAtomDipoleInitV(int n0,int l0, int l1) {

    if(l0 != l1 + 1 && l0 != l1 - 1) {
      string msg;
      msg = "|l0 - l1| != 1 in HAtomDipoleInitV";
      throw std::invalid_argument(msg);
    }

    typedef ExpBasis<F,1> Prim;
    typedef LinearComb<Prim> LC;
    
    LC psi_n = HAtomEigenState<F>(n0, l0);
    LC a = Op(OpDDr<Prim>(), psi_n);

    if( l0 > 0) {
      double coef;
      LC b;
      if(l0 < l1) 
	coef = - (l0 + 1);
      else 
        coef = l0;
      b = Op(OpCst<Prim>(coef), Op(OpRm<Prim>(-1), psi_n));
      a += b;
    }
    return a;    	
  }

  // ----------- eigen energy -------------------
  double HAtomEigenEnergy(int n) { return -1.0 / (2.0 * n * n); }
  
  
}
#endif



