#include <keys_values.hpp>
#include <l2func.hpp>
#include "factory.hpp"
#include "restrict.hpp"
#include "opt.hpp"
#include "opt_cbf.hpp"

namespace {
  using l2func::CSTO;
  using l2func::CGTO;
  using l2func::LinearComb;
  using l2func::HLikeAtom;
}
  
namespace opt_cbf_h {

  // ================ Utils ============================
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
  template<class Prim> 
  void checkBasis(const KeysValues& kv) {
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

  }

  // =============== Create Factory ====================
  string GetDrivType(const KeysValues& kv) {

    string di = kv.Get<string>("dipole");
    if(di == "length" || di == "velocity")
      return "STO";
    else if(di == "GTO")
      return "GTO";
    else {
      string msg; SUB_LOCATION(msg);
      throw InvalidDriv(msg, kv.Get<string>("channel"), di);
    }

  }
  IFactory* CreateFactory(const KeysValues& kv) {

    IFactory* res(NULL);

    int num_et  = kv.Count("opt_et_basis");
    int num_opt = kv.Count("opt_basis");
    string b_type = kv.Get<string>("basis_type");
    string d_type  = GetDrivType(kv);

    if( num_et == 0 && num_opt != 0) {
      if( b_type == "STO" && d_type == "STO")
	res = new FactoryMono<CSTO,CSTO>(kv);
      else if(b_type == "STO" && d_type == "GTO")
	res = new FactoryMono<CSTO,CGTO>(kv);
      else if(b_type == "GTO" && d_type == "STO")
	res = new FactoryMono<CGTO,CSTO>(kv);
      else if(b_type == "GTO" && d_type == "GTO")
	res = new FactoryMono<CGTO,CGTO>(kv);
      
    } else if( num_et != 0 && num_opt == 0) {
      if( b_type == "STO" && d_type == "STO")
	res = new FactoryEvenTemp<CSTO,CSTO>(kv);
      else if(b_type == "STO" && d_type == "GTO")
	res = new FactoryEvenTemp<CSTO,CGTO>(kv);
      else if(b_type == "GTO" && d_type == "STO")
	res = new FactoryEvenTemp<CGTO,CSTO>(kv);
      else if(b_type == "GTO" && d_type == "GTO")
	res = new FactoryEvenTemp<CGTO,CGTO>(kv);


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

  // =============== Interface ========================
  IFactory::IFactory(const KeysValues& _kv) : kv_(new KeysValues(_kv)) {}
  IFactory::~IFactory() {
    delete kv_;
  }

  template<class DrivPrim>
  void GetCustomDriv(const KeysValues& kv, LinearComb<DrivPrim>* driv) {

    int num = kv.Count("custom_driv");
    
    if(num == 0) {
      string msg; SUB_LOCATION(msg);
      msg += "\nkey named custom_driv not found.\n";
      throw runtime_error(msg);
    }

    typedef tuple<CD,int,CD> T;
    for(int i = 0; i < num; i++) {
      T d = kv.Get<T>("custom_driv");
      *driv += DrivPrim(get<0>(d), get<1>(d), get<2>(d));
    }
  }
  void GetChannel(const KeysValues& kv, HLikeAtom<CD>* h0, HLikeAtom<CD>* h1) {

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

    *h0 = HLikeAtom<CD>(n0, 1.0, l0);
    *h1 = HLikeAtom<CD>(-1, 1.0, l1); // -1 is dummy and meaning less

  }

  // static polymorphism for cSTO/cGTO is realized by using override
  void GetDriv(const KeysValues& kv, LinearComb<CSTO>* driv) {

    HLikeAtom<CD> h0;
    HLikeAtom<CD> h1;
    GetChannel(kv, &h0, &h1);

    string ch = kv.Get<string>("channel");
    string di = kv.Get<string>("dipole");
    int l1 = h1.l();
    
    if(di == "length") 
      *driv = h0.DipoleInitLength(l1);
    else if (di == "velocity") 
      *driv = h0.DipoleInitVelocity(l1);
    else if (di == "STO") {
      GetCustomDriv<CSTO>(kv, driv);
      string msg; SUB_LOCATION(msg); 
      msg += "\n mada junbi shiteinai \n";
      throw msg;
      *driv += 1.0 * l2func::CSTO(1.0, 2, 1.0);
    } else {
      string msg; SUB_LOCATION(msg); 
      throw InvalidDriv(msg, ch, di);
    }
  }
  void GetDriv(const KeysValues& kv, LinearComb<CGTO>* driv) {

    string ch = kv.Get<string>("channel");
    string di = kv.Get<string>("dipole");

    if (di == "GTO") {
      for(int i = 0; i < kv.Count("custom_driv"); i++) {
	typedef tuple<CD,int,CD> CIC;
	CIC cic = kv.Get<CIC>("custom_driv", i);
	CD  c = get<0>(cic);
	int n = get<1>(cic);
	CD  z = get<2>(cic);
	*driv += c * l2func::CGTO(1.0, n, z);
      }
    } else {
      string msg; SUB_LOCATION(msg); 
      throw InvalidDriv(msg, ch, di);      
    }
  }

  IOptTarget* IFactory::OptTarget() const {

    return this->optTarget();
    
  }
  int IFactory::BasisSize() const {
    
    VectorXcd zs;
    this->GetZs(&zs);

    return zs.rows();
  }

  // ============== Mono ================================
  template<class B, class D>
  FactoryMono<B,D>::FactoryMono(const KeysValues& kv) : IFactory(kv) {}
  template<class B, class D>
  FactoryMono<B,D>::~FactoryMono() {}

  template<class B, class D>
  vector<B>* FactoryMono<B,D>::BasisSet() const {

    checkBasis<B>(*kv_);

    vector<B>*  basis_set = new vector<B>();
    for(int i = 0; i < kv_->Count("opt_basis"); i++) {

      tuple<int, CD> n_z = kv_->Get<tuple<int, CD> >("opt_basis", i);
      B u(get<0>(n_z), get<1>(n_z), l2func::Normalized);
      basis_set->push_back(u);

    }
    
    return basis_set;

  }
  template<class B, class D>
  IOptTarget* FactoryMono<B,D>::optTarget() const {
    
    vector<B>* basis_set;
    l2func::HLikeAtom<CD> h_init;
    l2func::HLikeAtom<CD> h_final;
    l2func::LinearComb<D> driv;

    basis_set = this->BasisSet();
    GetChannel(*kv_, &h_init, &h_final);
    GetDriv(*kv_, &driv);
    double energy = kv_->Get<double>("energy");
    
    IOptTarget* ptr = new OptCBF<B,D>(*basis_set, h_final, driv, energy);
    return ptr;
  }
  template<class B, class D>
  IOptimizer<CD>* FactoryMono<B,D>::Optimizer() const {

    int max_iter = kv_->Get<int>("max_iter");
    double eps   = kv_->Get<double>("eps");
    IOptimizer<CD>* opt  = new OptimizerNewton<CD>(max_iter, eps);

    return opt;

  }
  template<class B, class D>
  void FactoryMono<B,D>::GetZs(VectorXcd* zs) const {

    vector<B>* us = this->BasisSet();
    *zs = VectorXcd::Zero(us->size());
    for(int i = 0; i < us->size(); i++)
	(*zs)[i] = (*us)[i].z();
    //    delete us;
    
  }

  // ============== EvenTempered ========================    
  template<class B, class D>
  FactoryEvenTemp<B,D>::FactoryEvenTemp(const KeysValues& kv) : IFactory(kv) {}
  template<class B, class D>
  FactoryEvenTemp<B,D>::~FactoryEvenTemp() {}
  template<class B, class D>
  vector<B>* FactoryEvenTemp<B,D>::BasisSet() const {

    checkBasis<B>(*kv_);
    
    vector<B>* basis_set = new vector<B>();
    int num_et = kv_->Count("opt_et_basis");
    
    for(int i = 0; i < num_et; i++) {
      int n, num_i;
      CD  z0, r;
      extractOptEtBasis2(*kv_, &n, &num_i, &z0, &r, i);      
      CD z = z0;
      for(int i = 0; i < num_i; i++) {
	B u(n, z, l2func::Normalized);
	basis_set->push_back(u);
	z *= r;
      }
    }
    return basis_set;
  }
  template<class B, class D>
  IOptTarget* FactoryEvenTemp<B,D>::optTarget() const {

    vector<B>* basis_set;
    l2func::HLikeAtom<CD> h_init;
    l2func::HLikeAtom<CD> h_final;
    l2func::LinearComb<D> driv;

    basis_set = this->BasisSet();
    GetChannel(*kv_, &h_init, &h_final);
    GetDriv(*kv_, &driv);
    double energy = kv_->Get<double>("energy");
    
    IOptTarget* ptr = new OptCBF<B,D>(*basis_set, h_final, driv, energy);
    delete basis_set;
    //    delete h_final; delete h_init; delete basis_set;
    return ptr;
  }
  template<class B, class D>
  IOptimizer<CD>* FactoryEvenTemp<B,D>::Optimizer() const {
    
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
  template<class B, class D>
  void FactoryEvenTemp<B,D>::GetZs(VectorXcd* zs) const {

    vector<B>* us = this->BasisSet();
    *zs = VectorXcd::Zero(us->size());
    for(int i = 0; i < us->size(); i++)
	(*zs)[i] = (*us)[i].z();
    delete us;
    
  }
}
