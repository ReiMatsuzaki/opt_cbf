#include "opt.hpp"
#include "restrict.hpp"
#include "from_kv.hpp"
#include <keys_values.hpp>
#include <l2func.hpp>
#include <macros.hpp>

using namespace l2func;

namespace opt_cbf_h {

  // --------- Basis Set -------------------
  void extractOptEtBasis(const KeysValues& kv, int* n, 
			 int* num, CD* z0, CD* r) {
    
    typedef tuple<int,int,CD,CD> IICC;
    IICC val = kv.Get<IICC>("opt_et_basis");
    *n   = get<0>(val);
    *num = get<1>(val);
    *z0  = get<2>(val);
    *r   = get<3>(val);
  }
  void errChk_BasisSet(int num_et, int num_opt) {

    if(num_et > 1) {
      string msg = "Multiple ET basis is not supported";
      throw runtime_error(msg);
    }
    if(num_et > 0 && num_opt > 0) {
      string msg = "# of opt_et_basis > 0 and # of opt_basis are not supported now";
      throw runtime_error(msg);
    }    

  }
  void err_BasisSet() {
    string msg = "Unsupported combination for basis";
      msg += " in Controller::setBasis\n";
      throw runtime_error(msg);
  }
  template<class Prim>
  void BuildBasisSet(const KeysValues& kv, 
		     vector<Prim>* basis_set) {
    
    int num_et = kv.Count("opt_et_basis");
    int num_opt= kv.Count("opt_basis");

    errChk_BasisSet(num_et, num_opt);
    
    vector<tuple<int, CD> > nz_list;
    int num;

    if(num_et == 1 && num_opt == 0) {
      int n;
      CD  z0, r;
      extractOptEtBasis(kv, &n, &num, &z0, &r);
      nz_list.resize(num);
      CD z = z0;
      for(int i = 0; i < n; i++) {
	nz_list[i] = make_tuple(n, z);
	z *= r;
      } 
    } else if (num_et == 0 && num_opt != 0) {
      num = kv.Count("opt_basis");
      nz_list.resize(num);
      typedef tuple<int, CD> I_CD;
      for(int i = 0; i < num; i++) 
	nz_list[i] = kv.Get<I_CD>("opt_basis", i);
    } else 
      err_BasisSet();

    for(int i = 0; i < num; i++) {
      int n = get<0>(nz_list[i]);
      CD  z = get<1>(nz_list[i]);
      Prim prim(n, z, Normalized);
      (*basis_set)[i] = prim;
    }

  }
  // ----------------------------------------

  /**
   * build driven term

  template<class Prim>
  void BuildHAtom(const KeysValues&, HAtomPI<Prim>*);
   */

  /**
   * build optimizer
   */
  void BuildOptimizer(const KeysValues& kv, 
		      IOptimizer<CD>* opt) {

    int max_iter = kv.Get<int>("max_iter");
    double eps   = kv.Get<double>("eps");

    // check number of basis 
    int num_et = kv.Count("opt_et_basis");
    int num_opt= kv.Count("opt_basis");

    if(num_et == 1 && num_opt == 0) {
      
      IRestriction<CD>* et;
      et = new EvenTemp<CD>();
      opt = new OptimizerRestricted<CD>
	(max_iter, eps, et);

    } else if(num_et == 0 && num_opt != 0) {

      opt = new OptimizerNewton<CD>(max_iter, eps, 0);

    } else {
      
      string msg = "Unsupported combination for basis";
      throw runtime_error(msg);

    }
  }

  

  // ---------- Explicit Instance --------------
  template void 
  BuildBasisSet<CSTO>(const KeysValues&, vector<CSTO>*);
  template void
  BuildBasisSet<CGTO>(const KeysValues&, vector<CGTO>*);
}
