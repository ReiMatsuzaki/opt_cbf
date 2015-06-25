#ifndef LCOMB_HPP
#define LCOMB_HPP

#include <vector>
#include <boost/function.hpp>
#include "prim.hpp"

namespace{
  using std::vector;
  using std::pair;
  using std::make_pair;
  using boost::function;
  using boost::bind;
}

namespace l2func {

  // =============== Linear Combination ===================
  template<class Prim>
  class LinearComb {

  public:
    // ----------- typedef -------------------
    typedef typename Prim::Field Field;
    typedef vector<pair<Field, Prim> > VFP;
    typedef typename VFP::const_iterator const_iterator;
    
  private:
    // ----------- Field Member ---------------
    VFP cf_list_; // c: coefficient, f:function
    
  public:
    // ---------- Constructor ---------------
    LinearComb();
    LinearComb(int n);
    LinearComb(const Prim& prim);
    
    // ---------- Accessor ------------------
    int size() const { return cf_list_.size(); }
    const_iterator begin() const { return cf_list_.begin(); }
    const_iterator end() const { return cf_list_.end(); }
    pair<Field, Prim>& operator [] (int i) {
      return cf_list_[i]; }
    const Prim& prim_i (int i) const {
      return cf_list_[i].second; }
    Field coef_i (int i) const { return cf_list_[i].first; }
    void operator += (pair<Field, Prim> cf);
    void operator += (const Prim& f);
    void operator += (const LinearComb<Prim>& o);
  };

  // =============== Functions ==================
  template<class Prim>
  pair<typename Prim::Field, Prim> operator *
  (typename Prim::Field c, const Prim& f) {
   
    return make_pair(c, f);
  }

  // --------------- Complex Inner Product ------
  template<class Prim1, class Prim2>
  typename Prim1::Field CIP(const LinearComb<Prim1>& a,
			    const LinearComb<Prim2>& b) {
    
    typedef typename LinearComb<Prim1>::const_iterator IT1;
    typedef typename LinearComb<Prim2>::const_iterator IT2;

    typename Prim1::Field acc(0);
    for(IT1 it_a = a.begin(), end_a = a.end(); it_a != end_a; ++it_a)
      for(IT2 it_b = b.begin(), end_b = b.end(); it_b != end_b; ++it_b) 
	acc += it_a->first * it_b->first * CIP(it_a->second, it_b->second);

    return acc;
  }
  template<class Prim, class L2Func>
  typename Prim::Field CIP(const LinearComb<Prim>& lc1,
			   const L2Func& o) {

    //    IsPrimitive<Prim>();

    typename Prim::Field acc(0);
    typedef typename LinearComb<Prim>::const_iterator IT;
    for(IT it = lc1.begin(), it_end = lc1.end(); it != it_end; ++it) {
      acc += it->first * CIP(it->second, o);
    }
    return acc;
  }
  template<class Prim, class L2Func>
  typename Prim::Field CIP( const L2Func& o,
			    const LinearComb<Prim>& lc1) {
    return CIP(lc1, o);
  }

  template<class Prim>
  LinearComb<Prim> Op(boost::function<Prim(const Prim&)> op,
		      const Prim& f) {
    typedef typename Prim::Field F;
    LinearComb<Prim> res;
    res += F(1) * op(f);
    return res;
  }
  template<class Prim>
  LinearComb<Prim> Op
  (boost::function<LinearComb<Prim>(const Prim&)> op,
   const Prim& f) {
    return op(f);
  }
  template<class Prim>
  LinearComb<Prim> Op(boost::function<Prim(const Prim&)> op,
		      const LinearComb<Prim>& f) {
    
    LinearComb<Prim> res;

    typedef typename LinearComb<Prim>::const_iterator IT;
    for(IT it = f.begin(), end_it = f.end();
	it != end_it; ++it) {

      typename Prim::Field ci = it->first;
      Prim op_fi = op(it->second);
      res += ci * op_fi;
      
    }

    return res;
  }
  template<class Prim>
  LinearComb<Prim> Op
  (boost::function<LinearComb<Prim>(const Prim&)> op,
		      const LinearComb<Prim>& f) {
    IsPrimitive<Prim>();
    LinearComb<Prim> g;

    typedef typename LinearComb<Prim>::const_iterator IT;
    for(IT it = f.begin(), it_end = f.end(); it != it_end; ++it) {

      typename Prim::Field ci = it->first;
      LinearComb<Prim> gi = op(it->second);

      for(IT it_i = gi.begin(); it_i != gi.end(); ++it_i) 
	g += (ci * it_i->first) * it_i->second;
    }

    return g;
  }

  template<class F, int m>
  LinearComb<ExpBasis<F,m> > OperateDDrForExp
  (const ExpBasis<F,m>& f) {
    
    typedef ExpBasis<F,m> Basis;

    F c = f.c();
    int n = f.n();
    F z = f.z();
    
    Basis a(c * F(n),     n-1, z);
    Basis b(c * (-z*F(m)), n+m-1, z);

    LinearComb<ExpBasis<F,m> > res;
    res += F(1) * a;
    res += F(1) * b;
    return(res);
  }
  template<class Prim>
  LinearComb<Prim> OperateDDr(const Prim& f) {
    IsPrimitive<Prim>();
    return OperateDDrForExp<typename Prim::Field,
			    Prim::exp_power>(f);
  }
  template<class Prim>
  boost::function<LinearComb<Prim>(const Prim&)> OpDDr() {
    return bind(OperateDDr<Prim>, _1);
  }

  template<class Prim>
  LinearComb<Prim> OperateDDr2 (const Prim& f) {
    LinearComb<Prim> df =  Op(OpDDr<Prim>(), f );
    LinearComb<Prim> ddf = Op(OpDDr<Prim>(), df);
    return ddf;
  }
  template<class Prim>
  boost::function<LinearComb<Prim>(const Prim&)> OpDDr2() {
    return bind(OperateDDr2<Prim>, _1);
  }

  template<class Prim>
  typename Prim::Field AtX(typename Prim::Field x,
			   const LinearComb<Prim>& f) {
    
    typename Prim::Field acc(0);

    typedef typename LinearComb<Prim>::const_iterator IT;
    for(IT it = f.begin(), end_it = f.end(); it != end_it; ++it) 
      acc += it->first * AtX(x, it->second);

    return acc;
  }

  // --------- derivative basis -------------
  // assuming a is normalized basis set.
  // above formula is from Sotsuron.
  // STO  
  // a(r) = N r^n Exp[-zr]
  //    N = 1/sqrt((2z)^{-2n-1} (2n)!) = z^{n+1/2} / sqrt(2n!)
  //   dN = (n+1/2)z^{n-1/2} / sqrt(2n!)
  //  ddN = (n-1/2)(n+1/2)z^{n-3/2} / sqrt(2n!)
  // dN/N = (n+1/2)/z
  //ddN/N = (n-1/2)(n+1/2)/(z^2)
  // GTO
  // a(r) = N r^n Exp[-zr^2]
  // dN/N = (1/4 + n/2) / z
  //ddN/N = (-3/4+n/2)(1/4+n/2)/(z^2)

  template<class F>
  LinearComb<ExpBasis<F,1> > D1Normalized(const ExpBasis<F, 1>& a) {
    LinearComb<ExpBasis<F, 1> > res;
    F   c = a.c(); // normalization term
    int n = a.n(); // principle number
    F   z = a.z(); // orbital expoent

    F   cp = (n + 0.5) / z * c;
    
    res += F(1) * ExpBasis<F, 1>(cp, n,   z);
    res += F(1) * ExpBasis<F, 1>(-c, n+1, z);

    return res;
  }
  template<class F>
  LinearComb<ExpBasis<F,2> > D1Normalized(const ExpBasis<F, 2>& a) {
    LinearComb<ExpBasis<F, 2> > res;
    F   c = a.c(); // normalization term
    int n = a.n(); // principle number
    F   z = a.z(); // orbital expoent

    F   cp = (n * 0.5 + 0.25) / z * c;
    
    res += 1.0 * ExpBasis<F, 2>(cp, n,   z);
    res += 1.0 * ExpBasis<F, 2>(-c, n+2, z);

    return res;
  }
  template<class F>
  LinearComb<ExpBasis<F,1> > D2Normalized(const ExpBasis<F, 1>& a) {
    LinearComb<ExpBasis<F, 1> > res;
    F   c = a.c(); // normalization term
    int n = a.n(); // principle number
    F   z = a.z(); // orbital expoent

    F   cp = (n + 0.5) / z * c;
    F  cpp = F(4 * n * n - 1) / (F(4) * z * z) * c;
    
    res += F(1) * ExpBasis<F, 1>(cpp,        n,   z);
    res += F(1) * ExpBasis<F, 1>(-F(2) * cp, n+1, z);
    res += F(1) * ExpBasis<F, 1>(c,          n+2, z);

    return res;
  }
  template<class F>
  LinearComb<ExpBasis<F,2> > D2Normalized(const ExpBasis<F, 2>& a) {
    LinearComb<ExpBasis<F, 2> > res;
    F   c = a.c(); // normalization term
    int n = a.n(); // principle number
    F   z = a.z(); // orbital expoent

    F   cp = (n * 0.5 + 0.25) / z * c;
    F  cpp = (-3.0/4.0 + n / 2.0) * (0.25 + n * 0.5) / (z*z) * c;    
    
    res += 1.0 * ExpBasis<F, 2>(cpp,       n,   z);
    res += 1.0 * ExpBasis<F, 2>(-F(2)*cp , n+2, z);
    res += 1.0 * ExpBasis<F, 2>(c,         n+4, z);

    return res;
  }    
  
}

#endif
