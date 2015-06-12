#ifndef OP_HPP
#define OP_HPP


/*
#include <vector>
#include <boost/functio.hpp>
#include "lcomb.hpp"

namespace {
  using std::vector;
  using std::pair;
  using std::make_pair;
  
  using boost::function;
}

namespace l2func {

  template<class Prim>
  class LinearOp {

  public:
    // ---------- typedef ----------------
    typedef function<LinearComb<Prim>(const Prim&)> OP;
    typedef Prim::Field Field;
    typedef vector<pair<Field, OP> >::const_iterator const_iterator;
    
  private:
    // ---------- Member Field ----------
    vector<pair<Field, OP> > coef_op_list_;

  public:
    // ------- Constructor --------------
    LinearOp() { IsPrimitive<Prim>(); }
    LinearOp(OP op) {
      coef_op_list_.resize(1);
      coef_op_list_[0] = make_pair(Field(1), op);
    }
    
    // ------- Getter -------------------
    int size() const { return coef_op_list_.size(); }
    const_iterator begin() const { return coef_op_list_.begin(); }
    const_iterator end()   const { return coef_op_list_.end(); }
    
    // ------- Setter ------------------
    void operator += (const pair<Field, Prim>& coef_op) {
      coef_op_list_.push_back(coef_op);
    }
    void operator += (const LinearOp& o) {
      
      for(const_iterator it = o.begin(); it != o.end(); ++it) 
	coef_op_list_.push_back(*it);
      
    }
    pair<Field, Prim>& operator [] (int i) { return coef_op_list_[i]; }

    // ------- Operation ---------------
    LinearComb<Prim> operator() (const Prim& a) {
      LinearComb<Prim> res;
      for(const_iterator it = coef_op_list_.begin(),
	    it_end = coef_op_list_.end(); it != it_end; ++it) {
	Field c = it->first;
	OP    f = it->second;
	res += c * f(a);
      }
      return res;
  };
}
*/

#endif
