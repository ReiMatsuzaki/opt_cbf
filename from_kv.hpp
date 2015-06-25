#ifndef FROM_KV_HPP
#define FROM_KV_HPP

/**
 * This file contains functions which calculate various objects
 * necessary in opt_cbf program from KeysValues.
 */

namespace opt_cbf_h {

  // --------- forward declaration ----------
  class KeysValues;
  template<class F> class HLikeAtom<F>;
  template<class Prim> class LinearComb<Prim>;

  /**
   * build basis set 
   */
  template<class Prim>
  void BasisSet(const KeysValues&, vector<Prim>*);
  
  /**
   * build driven term
   */
  template<class Prim>
  void BuildHAtom(const KeysValues&, HAtomPi<Prim>*);

  /**
   * build optimizer
   */
  void BuildOptimizer(const KeysValues&, IOptimizer*);
}

#endif
