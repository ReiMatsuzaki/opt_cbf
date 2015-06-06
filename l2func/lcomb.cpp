#include "lcomb.hpp"

namespace l2func {

  template<class Prim>
  LinearComb<Prim>::LinearComb() { IsPrimitive<Prim>();}
  template<class Prim>
  LinearComb<Prim>::LinearComb(int n):
    cf_list_(n) {}
  template<class Prim>
  LinearComb<Prim>::LinearComb(const Prim& prim):
    cf_list_(1) {
    cf_list_[0] = make_pair(1.0, prim);
  }

  template<class Prim>
  void LinearComb<Prim>::operator += (pair<Field, Prim> cf) {
    cf_list_.push_back(cf);
  }

  template<class Prim>
  void LinearComb<Prim>::operator += (const Prim& f) {
    cf_list_.push_back(make_pair(Field(1), f));
  }
  template<class Prim>
  void LinearComb<Prim>::operator += (const LinearComb<Prim>& o) {
    for(const_iterator it = o.begin(), it_end = o.end();
	it != it_end; ++it ) {
      cf_list_.push_back(*it);
    }
  }

  // explicit instance
  template class LinearComb<RSTO>;
  template class LinearComb<CSTO>;
  template class LinearComb<RGTO>;
  template class LinearComb<CGTO>;
}
