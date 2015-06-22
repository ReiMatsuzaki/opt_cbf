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
    virtual ~IRestriction();
    // from variables in normal space, set variables in rest.
    virtual void SetVars(const VecF& xs) = 0;
    // compute variables in normal space
    virtual VecF Xs() const = 0;
    // return the number of variables in restricted space
    virtual int size() const = 0;
    // compute Gradient of restricted variables
    virtual VecF Grad(const VecF&) const = 0;
    // compute Hessian of restricted variables
    virtual MatF Hess(const VecF&, const MatF&) const = 0;
    // shift vars in restricted space
    virtual void Shift(const VecF&) = 0;
  };

  template<class F>
  class EvenTemp : public IRestriction<F> {
    // ------- Typedef ----------------------
    typedef Matrix<F, Dynamic, 1>       VecF;
    typedef Matrix<F, Dynamic, Dynamic> MatF;
    // ------- Field ------------------------
    int num_;
    F x0_;
    F ratio_;
    
  public:
    // ------- Constructors -----------------
    EvenTemp();
    EvenTemp(int _num, F _r, F _x0);
    ~EvenTemp();
    // ------- Accessor ---------------------
    int size() const { return 2; }
    F ratio() const { return ratio_; }
    F x0()    const { return x0_; }
    // ------- Methods ----------------------
    void SetVars(const VecF&);
    VecF Xs() const;
    VecF Grad(const VecF&) const;
    MatF Hess(const VecF&, const MatF&) const;    
    void Shift(const VecF&);
  };
}


#endif
