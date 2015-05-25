#ifndef LCOMB_HPP
#define LCOMB_HPP

#include <vector>
#include <tr1/functional>

namespace{
  using std::vector;
  using std::pair;
  using std::make_pair;

  using std::tr1::function;
  using std::tr1::bind;
  using namespace std::tr1::placeholders;
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
    LinearComb() { IsPrimitive<Prim>();}
    LinearComb(int n):
      cf_list_(n) {}
    LinearComb(const Prim& prim):
      cf_list_(1) {
      cf_list_[0] = make_pair(1.0, prim);
    }
    // ---------- Accessor ------------------
    int size() const { return cf_list_.size(); }
    const_iterator begin() const { return cf_list_.begin(); }
    const_iterator end() const { return cf_list_.end(); }
    void operator += (pair<Field, Prim> cf) {
      cf_list_.push_back(cf);
    }
    void operator += (const Prim& f) {
      cf_list_.push_back(make_pair(Field(1), f));
    }
    void operator += (const LinearComb<Prim>& o) {
      for(const_iterator it = o.begin(), it_end = o.end();
	  it != it_end; ++it ) {
	cf_list_.push_back(*it);
      }
    }
    pair<Field, Prim>& operator [] (int i) { return cf_list_[i]; }
    const Prim& prim_i (int i) const { return cf_list_[i].second; }
    Field coef_i (int i) const { return cf_list_[i].first; }
  };

  // =============== Functions ==================
  template<class Prim>
  pair<typename Prim::Field, Prim> operator *(typename Prim::Field c,
					      const Prim& f) {
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

  // ---------------- operation ------------------
  /*  
  template<class Prim> class DDr1 {
  public:
    DDr1() {
      IsPrimitive<Prim>();
    }
    LinearComb<Prim> operator() (const Prim& prim) {
      
    }
  };
  template<class Prim> class OpRm {
  public:
    int m;
    OpRm(int _m) : m(_m) {}
    LinearComb<Prim> operator() (const Prim& a) {
      Prim prim(a.c(), a.n() + );
      return LinearComb<Prim>(
    }
  };
  */
  template<class Prim>
  LinearComb<Prim> Op(function<Prim(const Prim&)> op,
		      const Prim& f) {
    typedef typename Prim::Field F;
    LinearComb<Prim> res;
    res += F(1) * op(f);
    return res;
  }
  template<class Prim>
  LinearComb<Prim> Op(function<LinearComb<Prim>(const Prim&)> op,
		      const Prim& f) {
    return op(f);
  }
  template<class Prim>
  LinearComb<Prim> Op(function<Prim(const Prim&)> op,
		      const LinearComb<Prim>& f) {
    
    LinearComb<Prim> res;

    typedef typename LinearComb<Prim>::const_iterator IT;
    for(IT it = f.begin(), end_it = f.end(); it != end_it; ++it) {

      typename Prim::Field ci = it->first;
      Prim op_fi = op(it->second);
      res += ci * op_fi;
      
    }

    return res;
  }
  template<class Prim>
  LinearComb<Prim> Op(function<LinearComb<Prim>(const Prim&)> op,
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
  LinearComb<ExpBasis<F,m> > OperateDDrForExp(const ExpBasis<F,m>& f) {
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
    return OperateDDrForExp<typename Prim::Field, Prim::exp_power>(f);
  }
  template<class Prim>
  function<LinearComb<Prim>(const Prim&)> OpDDr() {
    return bind(OperateDDr<Prim>, _1);
  }

  template<class Prim>
  LinearComb<Prim> OperateDDr2 (const Prim& f) {
    LinearComb<Prim> df =  Op(OpDDr<Prim>(), f );
    LinearComb<Prim> ddf = Op(OpDDr<Prim>(), df);
    return ddf;
  }
  template<class Prim>
  function<LinearComb<Prim>(const Prim&)> OpDDr2() {
    return bind(OperateDDr2<Prim>, _1);
  }
}

#endif
