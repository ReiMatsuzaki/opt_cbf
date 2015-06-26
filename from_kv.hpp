#ifndef FROM_KV_HPP
#define FROM_KV_HPP

#include <complex>
#include <vector>

namespace {
  using std::vector;
  typedef std::complex<double> CD;
}

/**
 * This file contains functions which calculate various objects
 * necessary in opt_cbf program from KeysValues.
 */

// --------- forward declaration ----------
class KeysValues;
namespace l2func {
  template<class Prim> class LinearComb;
}
namespace opt_cbf_h {
  template<class F> class HAtomPI;
  template<class F> class IOptimizer;
}
namespace opt_cbf_h {
  
  /**
   * build basis set 
   */
  template<class Prim>
  void BuildBasisSet(const KeysValues&, vector<Prim>*);
  
  /**
   * build driven term
   */
  template<class Prim>
  void BuildHAtom(const KeysValues&, HAtomPI<Prim>**);

  /**
   * build optimizer
   */
  void BuildOptimizer(const KeysValues&, 
		      IOptimizer<CD>**);
}

#endif
