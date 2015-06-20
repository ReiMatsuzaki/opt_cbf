#ifndef RESTRICT_HPP
#define RESTRICT_HPP

#include <Eigen/Core>

namespace opt_cbf {

  template<class F>
  class IRestriction {
    typedef Matrix<F, Dynamic, 1>       VecF;
    typedef Matrix<F, Dynamic, Dynamic> MatF;
  public:
    virtual ~IRestriction();
    virtual void SetVars(const VecF& xs);
    virtual VecF Grad(const VecF&);
    virtual MatF Hess(const MatF&);
  };

  template<class F>
  class EvenTemp {
    // ------- Typedef ----------------------
    typedef Matrix<F, Dynamic, 1>       VecF;
    typedef Matrix<F, Dynamic, Dynamic> MatF;
    // ------- Field ------------------------
    F ratio_;
    F x0_;
    
  public:
    // ------- Methods ----------------------
    void SetVars(const VecF&);
    VecF Grad(const VecF&);
    MatF Hess(const MatF&);    
  };
}


#endif
