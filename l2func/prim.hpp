#ifndef STO_HPP
#define STO_HPP

#include <string>
#include <complex>
#include <tr1/functional>

namespace {
  using std::string;
  typedef std::complex<double> CD;
  using std::tr1::function;
  using std::tr1::bind;
  using namespace std::tr1::placeholders;
}

namespace l2func {

  // ==================== static ======================
  enum ENormalized {
    Normalized
  };
  template<typename T> struct raise_error;
  
  // ============= inner product ==========
  const CD operator*(const CD& a, int b);
  const CD operator*(int b, const CD& a);
  int int_pow(int base, unsigned int expo);
  template<class F> F STO_Int(F z, int n);
  template<class F> F GTO_Int(F z, int n);
       
  // ==================== STO or GTO =================
  // represent STO or GTO.
  // m==1 => STO, m==2 => GTO
  // value of m can be called exp_power;
  template<class F, int m>
  class ExpBasis {
    
  public:
    // ---------- typedef ---------------------------
    typedef F Field;
    enum EExpPower { exp_power=m };
    
    // ---------- Field Member ----------------------
  private:
    F c_;    // coefficient
    int n_;  // principle number
    F z_;    // orbital exponent
    //    bool is_zero_; // if true this basis is 0 in Hilbert space.
    
  public:
    // ----------- Constructors ---------------------
    ExpBasis();
    ExpBasis(int _n, F _z);
    ExpBasis(F _c, int _n, F _z);
    ExpBasis(int _n, F _z, ENormalized);
    template<class U>
    ExpBasis(const ExpBasis<U, m>& o):
      c_(o.c()), n_(o.n()), z_(o.z()) {}      
    // ----------- Accessors ------------------------
    F c() const { return c_; }
    int n() const { return n_; }
    F z() const { return z_; }
    void set_z(F z) { z_ = z; }
  };

  // =========== typedef =========================
  
  typedef ExpBasis<double, 1> RSTO;
  typedef ExpBasis<double, 2> RGTO;
  typedef ExpBasis<CD, 1> CSTO;
  typedef ExpBasis<CD, 2> CGTO;
  
  template<class F>
  F CIP(const ExpBasis<F, 1>& a, const ExpBasis<F, 1>& b);
  template<class F>
  F CIP(const ExpBasis<F, 1>& a, const ExpBasis<F, 2>& b);
  template<class F>  
  F CIP(const ExpBasis<F, 2>& a, const ExpBasis<F, 1>& b);
  template<class F>
  F CIP(const ExpBasis<F, 2>& a, const ExpBasis<F, 2>& b);

  // =========== raise error if not primitive =========
  template<class Prim> struct IsPrimitive;
  template<> struct IsPrimitive<RSTO> {};
  template<> struct IsPrimitive<CSTO> {};
  template<> struct IsPrimitive<RGTO> {};
  template<> struct IsPrimitive<CGTO> {};

  // ============= operation =================
  template<class Prim>
  typename Prim::Field AtX(typename Prim::Field x, const Prim& f) {
    int m = Prim::exp_power;
    return f.c() * pow(x, f.n()) * exp(-f.z() * pow(x, m));
  }
  template<int num, class Prim>
  Prim DBasis(typename Prim::Field c, int n,
	      typename Prim::Field z) {
    
    IsPrimitive<Prim>();
    int m = Prim::exp_power;
    return Prim(pow(-1.0, num) * c, n + m * num, z);
  }
  template<class Prim>
  Prim OperateRm( int m, const Prim& f) {
    return Prim(f.c(), f.n() + m, f.z());
  }
  template<class Prim>
  Prim OperateCst(typename Prim::Field c, const Prim& f) {
    return Prim(f.c() * c, f.n(), f.z());
  }
  template<class Prim>
  function<Prim(const Prim&)> OpRm(int m) {
    return bind(OperateRm<Prim>, m, _1);
  }
  template<class Prim>
  function<Prim(const Prim&)> OpCst(typename Prim::Field c) {
    return bind(OperateCst<Prim>, c, _1);
  }  
}

#endif
