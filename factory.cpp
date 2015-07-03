#include <keys_values.hpp>
#include <l2func.hpp>
#include "factory.hpp"
#include "driv.hpp"
#include "restrict.hpp"
#include "opt.hpp"

namespace opt_cbf_h {

  // =============== Create Factory ====================
  IFactory* CreateFactory(const KeysValues& kv) {

    IFactory* res(NULL);

    int num_et  = kv.Count("opt_et_basis");
    int num_opt = kv.Count("opt_basis");

    if( num_et == 0 && num_opt != 0) {
      
      res = NULL;
      
    } else if( num_et != 0 && num_opt == 0) {

      res = new FactoryEvenTemp(kv);


    } else if( num_et == 0 && num_opt == 0) {

      string msg; SUB_LOCATION(msg); 
      msg += "\nnumber of basis is zero\n";
      throw runtime_error(msg);      

    } else {

      std::string msg; SUB_LOCATION(msg);
      msg += "\n# of opt_et_basis = 0 or # of opt_basis = 0\n";
      throw runtime_error(msg);

    }

    return res;    
  }

  // ============== Exceptions =========================
  InvalidBasis::InvalidBasis(string msg, string basis_type) :
    runtime_error(msg), basis_type_(basis_type) {}
  InvalidBasis::~InvalidBasis() throw() {}
  const char* InvalidBasis::what() const throw() {
    string msg = runtime_error::what();
    msg += "\ninvalid basis type.\n";
    msg += "basis_type : ";
    msg += basis_type_;
    return msg.c_str();
  }
  InvalidDriv::InvalidDriv(string msg, string ch, string di) : 
    runtime_error(msg), ch_(ch), di_(di){}
  InvalidDriv::~InvalidDriv() throw() {}
  const char* InvalidDriv::what() const throw() {
    string msg(runtime_error::what());
    msg += "\ninvalid dipole or channel\n";
    msg += "\nchannel: ";
    msg += ch_;
    msg += "\n";
    msg += "dipole: ";
    msg += di_;
    return msg.c_str();
  }

  // =============== Abstract Interface ================
  IFactory::IFactory(const KeysValues& _kv) : kv_(new KeysValues(_kv)) {}
  IFactory::~IFactory() {
    delete kv_;
  }

  // ============== EvenTempered ========================  
  template<class Prim> string basisName();
  template<> string basisName<l2func::CSTO>() { return "STO"; }
  template<> string basisName<l2func::CGTO>() { return "GTO"; }
  void extractOptEtBasis2(const KeysValues& kv, int* n, int* num, CD* z0, 
			 CD* r, int idx) {
			 
    typedef tuple<int,int,CD,CD> IICC;
    IICC val = kv.Get<IICC>("opt_et_basis", idx);
    *n   = get<0>(val);
    *num = get<1>(val);
    *z0  = get<2>(val);
    *r   = get<3>(val);

  }
  
  FactoryEvenTemp::FactoryEvenTemp(const KeysValues& kv) : IFactory(kv) {}
  FactoryEvenTemp::~FactoryEvenTemp() {}

  template<class Prim> vector<Prim>* basisSet(const KeysValues& kv) {

    string basis_type;
    try {
      basis_type = kv.Get<string>("basis_type");
    } catch(exception& e) {
      throw runtime_error("\nbasis_type does not found\n");
    }

    if(kv.Get<string>("basis_type") != basisName<Prim>()) {
      string msg; SUB_LOCATION(msg); 
      throw InvalidBasis(msg, basisName<Prim>());
    }
    
    int num_et = kv.Count("opt_et_basis");
    vector<Prim>* basis_set = new vector<Prim>();
    for(int i = 0; i < num_et; i++) {
      int n, num_i;
      CD  z0, r;
      extractOptEtBasis2(kv, &n, &num_i, &z0, &r, i);      
      CD z = z0;
      for(int i = 0; i < num_i; i++) {
	Prim u(n, z, l2func::Normalized);
	basis_set->push_back(u);
	z *= r;
      }
    }
    return basis_set;
}
  vector<l2func::CSTO>* FactoryEvenTemp::STOSet() const {
    return basisSet<l2func::CSTO>(*kv_);
  }
  vector<l2func::CGTO>* FactoryEvenTemp::GTOSet() const {
    return basisSet<l2func::CGTO>(*kv_);
  }

  template<class Prim> HAtomPI<Prim>* hAtomPI(const KeysValues& kv) {

    // check basis type
    string basis_type = kv.Get<string>("basis_type");
    if(basis_type != basisName<Prim>()) {
      string msg; SUB_LOCATION(msg); 
      throw InvalidBasis(msg, basis_type);
    }

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
    } else {
      l0 = -1; l1 = -1; n0 = -1;
      string msg; SUB_LOCATION(msg); 
      throw InvalidDriv(msg, ch, di);
    }

    // driven term
    l2func::HLikeAtom<CD> hatom(n0, 1.0, l0);
    l2func::LinearComb<l2func::CSTO> mu_phi;
    if(di == "length")
      mu_phi = hatom.DipoleInitLength(l1);
    else if (di == "velocity")
      mu_phi = hatom.DipoleInitVelocity(l1);
    else  {
      string msg; SUB_LOCATION(msg); 
      throw InvalidDriv(msg, ch, di);
    }

    // energy 
    double ene = kv.Get<double>("energy");

    // create HAtomPI object
    HAtomPI<Prim>* h_pi = new HAtomPI<Prim>(l1, 1.0, ene, mu_phi);
    return h_pi;
  }
  HAtomPI<l2func::CSTO>* FactoryEvenTemp::HAtomPiSTO() const {
    return hAtomPI<l2func::CSTO>(*kv_);
  }
  HAtomPI<l2func::CGTO>* FactoryEvenTemp::HAtomPiGTO() const {
    return hAtomPI<l2func::CGTO>(*kv_);
  }

  IOptimizer<CD>* FactoryEvenTemp::Optimizer() const {
    
    int max_iter = kv_->Get<int>("max_iter");
    double eps   = kv_->Get<double>("eps");

    vector<int> num_list;
    for(int idx = 0; idx < kv_->Count("opt_et_basis"); idx++) {
      int n, num_i;
      CD  z0, r;
      extractOptEtBasis2(*kv_, &n, &num_i, &z0, &r, idx);
      num_list.push_back(num_i);
    }

    IRestriction<CD>* et = new MultiEvenTemp<CD>(num_list);
    IOptimizer<CD>* opt  = new OptimizerRestricted<CD>(max_iter, eps, et);

    return opt;
  }
}
