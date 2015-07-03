#ifndef RESTRICT_HPP
#define RESTRICT_HPP

#include <Eigen/Core>
#include <vector>
#include <boost/tuple/tuple.hpp>

namespace {
  using namespace Eigen;
  using std::vector;
  using boost::tuple;
  using boost::make_tuple;
}

namespace opt_cbf_h {

  // ============= Interface ================
  template<class F>
  class IRestriction {
  public:
    // --------- type -------------------
    typedef Matrix<F, Dynamic, 1>       VecF;
    typedef Matrix<F, Dynamic, Dynamic> MatF;

  public:
    // ---------- methods ---------------
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

  // ============== No restriction ==========
  template<class F>
  class NoRestriction : public IRestriction<F> {
  public:
    // ------- Typedef ----------------------
    typedef Matrix<F, Dynamic, 1>       VecF;
    typedef Matrix<F, Dynamic, Dynamic> MatF;

  private:
    // ------- Field ------------------------
    VecF xs_;

  public:
    // ------- Constructors -----------------
    // NoRestriction(VecF _xs) : xs_(_xs) {}
    NoRestriction() {}
    ~NoRestriction() {}
    // ------- Accessor --------------------
    int size() const { return xs_.rows(); }
    // ------- Methods ---------------------
    void SetVars(const VecF& xs) {
      xs_ = VecF::Zero(xs.rows());
      xs_ = xs; }
    VecF Xs() const { return xs_; }
    VecF Grad(const VecF& g) const { return g; }
    MatF Hess(const VecF&, const MatF& h) const { return h; }
    void Shift(const VecF& dx) { xs_ += dx; }
  };

  // ============== EvenTempered ============
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
    // EvenTemp(int _num, F _r, F _x0);
    ~EvenTemp();
    // ------- Accessor ---------------------
    int size() const { return num_; }
    F ratio() const { return ratio_; }
    F x0()    const { return x0_; }
    // ------- Methods ----------------------
    void SetVars(const VecF&);
    VecF Xs() const;
    VecF Grad(const VecF&) const;
    MatF Hess(const VecF&, const MatF&) const;    
    void Shift(const VecF&);
  };

  // ============== Multiple ET =============
  template<class F>
  class MultiEvenTemp : public IRestriction<F> {
  public:
    // -------- typedef -----------
    typedef typename IRestriction<F>::VecF VecF;
    typedef typename IRestriction<F>::MatF MatF;
    typedef tuple<int,F,F> IFF;
    typedef typename vector<IFF>::iterator IT;
    typedef typename vector<IFF>::const_iterator CIT;
  
  private:
    // -------- Field Member ------
    // num_list_ = [(2,a,r), (3,b,s), (3,c,t)] means that
    // x0=a, x1=ar
    // x2=b, x3=bs, x4=bss
    // x5=c, x6=ct, x7=ctt
    vector<IFF> num_x0_r_list_;

  public:
    // -------- Constructors ------
    MultiEvenTemp(const vector<int>& index_list);
    ~MultiEvenTemp();

    // ------- Accessor -----------
    IFF num_x0_r(int i) const;

    // -------- Methods -----------
    void SetVars(const VecF& xs);
    VecF Xs() const;
    int size() const;
    VecF Grad(const VecF&) const;
    MatF Hess(const VecF&, const MatF&) const;
    void Shift(const VecF&);
    
  };
}


#endif
