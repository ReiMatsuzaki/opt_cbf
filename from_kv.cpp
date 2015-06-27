#include "opt.hpp"
#include "restrict.hpp"
#include "from_kv.hpp"
#include "opt.hpp"
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
      for(int i = 0; i < num; i++) {
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

    basis_set->resize(num);
    for(int i = 0; i < num; i++) {
      int n = get<0>(nz_list[i]);
      CD  z = get<1>(nz_list[i]);
      Prim prim(n, z, Normalized);
      (*basis_set)[i] = prim;
    }

  }
  // ----------------------------------------

  // --------- HAtom Photoionization -----------
  void err_BuildHAtomPI1(string ch) {
    string msg;
    msg = "unsupported channel\n";
    msg+= "channel: ";
    msg+= ch;
    throw runtime_error(msg);
  }
  void err_BuildHAtomPI2(string di) {
    string msg;
    msg = "unsupported dipole operator.\n";
    msg += "dipole : ";
    msg += di;
    throw runtime_error(msg);
  }
  void err_BuildHAtomPI3(string basis_type) {
    string msg;
    msg =  "invalid basis type\n";
    msg+= "basis_type: ";
    msg+= basis_type;
    throw runtime_error(msg);
  }
  template<class Prim>
  void BuildHAtomPI1(const KeysValues& kv, 
		     HAtomPI<Prim>** hatom) {

    // Hydrogen atom
    string ch = kv.Get<string>("channel");
    string di = kv.Get<string>("dipole");
    int l0, l1, n0;
    if(ch == "1s->kp") {
      l0 = 0; l1 = 1; n0 = 1;
    } else if(ch == "2p->ks") {
      l0 = 1; l1 = 0; n0 = 2;
    } else if(ch == "2p->kd") {
      l0 = 1; l1 = 2; n0 = 2;
    } else if(ch == "3d->kp") {
      l0 = 2; l1 = 1; n0 = 3;
    } else if(ch == "3d->kf") {
      l0 = 2; l1 = 3; n0 = 3;
    } else 
      err_BuildHAtomPI1(ch);
      
    double energy = kv.Get<double>("energy");
    HLikeAtom<CD> hatom(n0, 1.0, l0);
    LinearComb<CSTO> mu_phi;

    if(di == "length")
      mu_phi = hatom.DipoleInitLength(l1);
    else if (di == "velocity")
      mu_phi = hatom.DipoleInitVelocity(l1);
    else 
      err_BuildHAtomPI2(di);
     
    string basis_type = 
      kv.Get<string>("basis_type");
    IOptTarget* ptr;

    if(basis_type == "STO") {
      vector<CSTO> basis_set;
      this->setBasis<CSTO>(&basis_set);
      *hatom = new HAtomPI<CSTO>(l1, 1.0, energy, mu_phi);
    } else if (basis_type == "GTO") {
      vector<CGTO> basis_set;
      this->setBasis<CGTO>(&basis_set);
      *hatom = new HAtomPI<CGTO>(l1, 1.0, energy, mu_phi);
    } else 
      err_BuildHAtomPI3();

  }
  // ---------------------------------------------

  // ----------- Optimizer --------------------
  void BuildOptimizer(const KeysValues& kv, 
		      IOptimizer<CD>** opt) {

    int max_iter = kv.Get<int>("max_iter");
    double eps   = kv.Get<double>("eps");

    // check number of basis 
    int num_et = kv.Count("opt_et_basis");
    int num_opt= kv.Count("opt_basis");

    if(num_et == 1 && num_opt == 0) {
      
      IRestriction<CD>* et;
      et = new EvenTemp<CD>();
      *opt = new OptimizerRestricted<CD>
	(max_iter, eps, et);

    } else if(num_et == 0 && num_opt != 0) {

      *opt = new OptimizerNewton<CD>(max_iter, eps, 0);

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
  template void
  BuildHAtomPI<CSTO>(const KeysValues&, HAtomPI<CSTO>**);
  template void
  BuildHAtomPI<CGTO>(const KeysValues&, HAtomPI<CGTO>**);
}
