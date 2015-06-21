#ifndef RESTRICT_HPP
#define RESTRICT_HPP

#include <Eigen/Core>

namespace {
  using namespace Eigen;
}

namespace opt_cbf_h {

  template<class F>
  class IRestriction {
    typedef Matrix<F, Dynamic, 1>       VecF;
    typedef Matrix<F, Dynamic, Dynamic> MatF;
  public:
    virtual ~IRestriction() = 0;
    virtual void SetVars(const VecF& xs) = 0;
    VecF Xs() const = 0;
    virtual VecF Grad(const VecF&) const = 0;
    virtual MatF Hess(const VecF&, const MatF&) const = 0;
  };

  template<class F>
  class EvenTemp {
    // ------- Typedef ----------------------
    typedef Matrix<F, Dynamic, 1>       VecF;
    typedef Matrix<F, Dynamic, Dynamic> MatF;
    // ------- Field ------------------------
    int num_;
    F ratio_;
    F x0_;
    
  public:
    // ------- Constructors -----------------
    EvenTemp();
    EvenTemp(int _num, F _r, F _x0);
    ~EvenTemp() {}
    // ------- Accessor ---------------------
    F ratio() const { return ratio_; }
    F x0()    const { return x0_; }
    // ------- Methods ----------------------
    void SetVars(const VecF&);
    VecF Xs() const;
    VecF Grad(const VecF&) const;
    MatF Hess(const VecF&, const MatF&) const;    
  };
}


#endif
