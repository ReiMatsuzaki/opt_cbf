#include <iostream>
#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include <keys_values.hpp>
#include <timer.hpp>
#include <l2func.hpp>
#include "controller.hpp"
#include "opt_cbf.hpp"
#include "driv.hpp"
#include "restrict.hpp"
#include "opt.hpp"

namespace {
  using std::cout;
  using std::endl;
  using std::ifstream;
  typedef std::complex<double> CD;
  using boost::shared_ptr;
  using namespace l2func;
}
namespace opt_cbf_h {

  class OptCBFController::Impl {
  public:
    // -------- type -----------------
    typedef tuple<int, CD> I_CD;
    typedef tuple<bool, string> BS;

    // -------- Member Field ---------
    KeysValues keys_values_;
    shared_ptr<IOptTarget> opt_target_;
    shared_ptr<IOptimizer<CD> > optimizer_;
    Timer      timer_;
    VectorXcd  zs_;
    OptRes<CD> opt_res_;

    // -------- Constructor ---------
    Impl() : keys_values_(":", " ") {}

    // ------- Read support ---------
    void convertWriteOption(const string& key) {

      BS default_val = make_tuple(false, "");
      keys_values_.SetIfNull<BS>(key, default_val);
      keys_values_.Check<bool, string>(NumberIs(1), key);

    }
    void convertData() {
      
      keys_values_.Check<string>(NumberIs(1), "channel");
      keys_values_.Check<string>(NumberIs(1), "dipole");
      keys_values_.Check<double>(NumberIs(1), "energy");

      keys_values_.Check<string>(NumberIs(1), "basis_type");
      keys_values_.Check<int,CD>(AnyNumber(), "opt_basis");
      keys_values_.Check<int,int,CD,CD>(AnyNumber(),
					"opt_et_basis"); 
					
      keys_values_.SetIfNull<int>("max_iter", 100);
      keys_values_.Check<int>(NumberIs(1), "max_iter");
      keys_values_.SetIfNull<double>("eps", 0.0000001);
      keys_values_.Check<double>(NumberIs(1), "eps");

      typedef tuple<double, double> DD;
      DD grid = make_tuple(10.0, 0.1);
      keys_values_.SetIfNull<DD>("grid", grid);
      keys_values_.Check<double, double>(NumberIs(1), "grid");

      this->convertWriteOption("write_psi");
      this->convertWriteOption("write_hess");
      this->convertWriteOption("write_grad");
    }
    template<class Prim>
    void setBasis(vector<Prim>* basis_set) {

      // check number of basis 
      int num_et = keys_values_.Count("opt_et_basis");
      int num_opt= keys_values_.Count("opt_basis");

      // error 
      if(num_et > 1) {
	string msg = "Multiple ET basis is not supported";
	throw runtime_error(msg);
      }

      if(num_et > 0 && num_opt > 0) {
	string msg = "# of opt_et_basis > 0 and # of opt_basis are not supported now";
	throw runtime_error(msg);
      }

      if(num_et == 1 && num_opt == 0) {
	typedef tuple<int,int,CD,CD> IICC;
	IICC val = keys_values_.Get<IICC>("opt_et_basis");
	int n   = get<0>(val);
	int num = get<1>(val);
	CD  x0  = get<2>(val);
	CD  r   = get<3>(val);
	zs_ = VectorXcd::Zero(num);

	CD z = x0;
	basis_set->resize(num);
	for(int i = 0; i < n; i++) {
	  z *= r;
	  Prim prim(n, z, Normalized);
	  (*basis_set)[i] = prim;
	  zs_(i) = z;
	}
      } else if(num_et == 0 && num_opt != 0) {

	int num = keys_values_.Count("opt_basis");
	basis_set->resize(num);
	zs_ = VectorXcd::Zero(num);

	for(int i = 0; i < num; i++) {
	
	  I_CD n_z = keys_values_.Get<I_CD>("opt_basis", i);
	  int n = get<0>(n_z);
	  CD  z = get<1>(n_z);
	  Prim prim(n, z, Normalized);
	  (*basis_set)[i] = prim;
	  zs_(i) = z;
	}
      } else {
	
	string msg = "Unsupported combination for basis";
	msg += " in Controller::setBasis\n";
	throw runtime_error(msg);
	
      }

    }
    void setOptTarget() {

      // Hydrogen atom
      string ch = keys_values_.Get<string>("channel");
      string di = keys_values_.Get<string>("dipole");
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
	string msg;
	msg = "unsupported channel\n";
	msg+= "channel: ";
	msg+= ch;
	throw runtime_error(msg);
      }
      double energy = keys_values_.Get<double>("energy");
      HLikeAtom<CD> hatom(n0, 1.0, l0);
      LinearComb<CSTO> mu_phi;

      if(di == "length")
	mu_phi = hatom.DipoleInitLength(l1);
      else if (di == "velocity")
	mu_phi = hatom.DipoleInitVelocity(l1);
      else {
	string msg;
	msg = "unsupported dipole operator.\n";
	msg += "dipole : ";
	msg += di;
	throw runtime_error(msg);
      }
     
      string basis_type = 
	keys_values_.Get<string>("basis_type");
      IOptTarget* ptr;

      if(basis_type == "STO") {
	vector<CSTO> basis_set;
	this->setBasis<CSTO>(&basis_set);
	HAtomPI<CSTO>* h_pi = 
	  new HAtomPI<CSTO>(l1, 1.0, energy, mu_phi);
	ptr = new OptCBF<CSTO>(basis_set, h_pi);
      } else if (basis_type == "GTO") {
	vector<CGTO> basis_set;
	this->setBasis<CGTO>(&basis_set);
	HAtomPI<CGTO>* h_pi = 
	  new HAtomPI<CGTO>(l1, 1.0, energy, mu_phi);
	ptr = new OptCBF<CGTO>(basis_set, h_pi);
      } else {
	string msg;
	msg =  "invalid basis type\n";
	msg+= "basis_type: ";
	msg+= basis_type;
	throw runtime_error(msg);
      }
      opt_target_ = shared_ptr<IOptTarget>(ptr);

    }
    void setOptimizer() {

      int max_iter = keys_values_.Get<int>("max_iter");
      double eps   = keys_values_.Get<double>("eps");
      IOptimizer<CD>* opt;

      // check number of basis 
      int num_et = keys_values_.Count("opt_et_basis");
      int num_opt= keys_values_.Count("opt_basis");

      if(num_et == 1 && num_opt == 0) {

	IRestriction<CD>* et;
	et = new EvenTemp<CD>();
	opt = new OptimizerRestricted<CD>(max_iter, eps, et);

      } else if(num_et == 0 && num_opt != 0) {

	opt = new OptimizerNewton<CD>(max_iter, eps, 0);

      } else {

	string msg = "Unsupported combination for basis";
	throw runtime_error(msg);

      }

      optimizer_ = shared_ptr<IOptimizer<CD> >(opt);

    }

    // -------- Write support --------
    void writeOptionalFiles() {

      // psi
      typedef tuple<double, double> DD;
      DD rr = keys_values_.Get<DD>("grid");
      double rmax = get<0>(rr);
      double dr   = get<1>(rr);
      BS bool_file1 = keys_values_.Get<BS>("write_psi");
      if(get<0>(bool_file1)) {
	string fn = get<1>(bool_file1);
	opt_target_->WritePsi(fn, rmax, dr);
      }
      
      // grad
      BS bool_file_grad = keys_values_.Get<BS>("write_grad");
      if(get<0>(bool_file_grad)) {
	string fn = get<1>(bool_file_grad);
	ofstream ofs(fn.c_str());
	if(ofs.fail()) {
	  string msg = "Failed to open file.\n";
	  msg += "File name : ";
	  msg += fn;
	  throw runtime_error(msg);
	}
	ofs << opt_res_.grad << endl;
      }
      
      // hess
      BS bool_file2 = keys_values_.Get<BS>("write_hess");
      if(get<0>(bool_file2)) {
	
	ofstream ofs(get<1>(bool_file2).c_str());
	if(ofs.fail()) {
	  string msg = "Failed to open file.\n";
	  msg += "File name : ";
	  msg += get<1>(bool_file2);
	  throw runtime_error(msg);
	}
	ofs << opt_res_.hess << endl;
      }
    }
    void writeBasis(ofstream& ofs) {
      
      int num_et = keys_values_.Count("opt_et_basis");
      int num_opt= keys_values_.Count("opt_basis");

      if(num_et == 0) {
	for(int i = 0; i < opt_res_.z.rows(); i++) {
	  ofs << "opt_basis: ";
	  ofs << get<0>(keys_values_.Get<I_CD>
			("opt_basis", i));
	  ofs << " " << opt_res_.z(i) << endl;
	}	
      }
      
      if(num_et == 1 && num_opt == 0) {
	typedef tuple<int,int,CD,CD> IICC;
	IICC val = keys_values_.Get<IICC>("opt_et_basis");
	int n   = get<0>(val);
	int num = get<1>(val);
	CD  x0  = get<2>(val);
	CD  r   = get<3>(val);
	CD z = x0;
	for(int i = 0; i < num; i++) {
	  ofs << "opt_basis: ";
	  ofs << n << " " << z << endl;
	  z *= r;
	}
      }
    }
    
    // -------- Method --------------
    void Read(const char* filename) {

      timer_.Start("read");
      
      ifstream ifs(filename);

      if(ifs.fail()) {
	string msg = "Failed to open file.\n";
	msg += "File name: ";
	msg += filename;
	throw runtime_error(msg);
      }

      keys_values_.Read(ifs);
      
      try {
	this->convertData();
      } catch (exception& e) {
	string msg = "failed to convertData. ";
	msg += "Error message is: \n";
	msg += e.what();
	throw runtime_error(msg);
      }
      
      this->setOptTarget();
      this->setOptimizer();

      timer_.End("read");
    }
    void Compute() {

      timer_.Start("calc");

      opt_res_ = optimizer_->Optimize
	(bind(&IOptTarget::Compute, opt_target_,
	      _1, _2, _3, _4),
	 zs_);
      
      timer_.End("calc");

    }
    void Write(const char* filename) {
      
      ofstream ofs(filename);

      if(ofs.fail()) {
	string msg = "Failed to open file.\n";
	msg += "File name : ";
	msg += filename;
	throw runtime_error(msg);
      }
      
      ofs << "convergence: ";
      ofs << (opt_res_.convergence ? "true" : "false") << endl;
      ofs << "iter_num: " << opt_res_.iter_num << endl;

      this->writeBasis(ofs);
      
      VectorXcd cs = opt_target_->GetCoefs();
      for(int i = 0; i < cs.rows(); i++) {
	ofs << "coef: ";
	ofs << cs(i) << endl;
      }
      ofs << "alpha: " << opt_res_.value << endl;
      ofs << "read_time: " << timer_.GetTime("read") << endl;
      ofs << "calc_time: " << timer_.GetTime("calc") << endl;

      this->writeOptionalFiles();
    }
  };
    
  // =========== interface =============
  OptCBFController::OptCBFController() : impl_(new Impl()) {}
  OptCBFController::~OptCBFController() { 
    delete impl_; }
  void OptCBFController::Read(const char* fn) { 
    impl_->Read(fn); }
  void OptCBFController::Compute() { impl_->Compute(); }
  void OptCBFController::Write(const char* fn) { 
    impl_->Write(fn); }
}
