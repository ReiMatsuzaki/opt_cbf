#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
#include <Eigen/Core>
#include <keys_values.hpp>
#include <timer.hpp>
#include <l2func.hpp>
#include "controller.hpp"
#include "opt_cbf.hpp"
#include "driv.hpp"
#include "restrict.hpp"
#include "opt.hpp"
#include "factory.hpp"

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
    boost::shared_ptr<IFactory> factory_;
    Timer      timer_;
    int debug_lvl_;
    boost::shared_ptr<IOptTarget> opt_target_;
    boost::shared_ptr<IOptimizer<CD> > optimizer_;
    VectorXcd  zs_;
    OptRes<CD> opt_res_;

    // -------- Constructor ---------
    explicit Impl(int _debug_lvl) : 
      keys_values_(":", " "), 
      timer_(),
      debug_lvl_(_debug_lvl)
    {}

    // ------- General support ------
    void PrintInDebug(string msg) {
      if(debug_lvl_ > 0)
	cout << msg << endl;
    }
    
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
    void setOptTarget() {

      IOptTarget* ptr(NULL);
      ptr = factory_->OptTarget();
      opt_target_ = boost::shared_ptr<IOptTarget>(ptr);
      factory_->GetZs(&zs_);

    }
    void setOptimizer() {
      
      IOptimizer<CD>* ptr(NULL);
      ptr = factory_->Optimizer();
      optimizer_ = boost::shared_ptr<IOptimizer<CD> >(ptr);

    }

    // ------- Compute support --------
    void err_OptimizerNull() {
      string msg;
      msg = "Optimizer is null\n";
      throw runtime_error(msg);
    }
    void err_OptTargetNull() {
      string msg;
      msg = "OptTarget is null\n";
      throw runtime_error(msg);
    }
    void err_zsNull() {
      string msg;
      msg = "zs is null\n";
      throw runtime_error(msg);
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
      
      if(num_opt == 0) {
	typedef tuple<int,int,CD,CD> IICC;
	int num_et = keys_values_.Count("opt_et_basis");
	int acc_num(0);
	for(int i = 0; i < num_et; i++) {
	  IICC val = keys_values_.Get<IICC>("opt_et_basis", i);
	  int n   = get<0>(val);
	  int num = get<1>(val);
	  CD  x0  = opt_res_.z(acc_num);
	  CD  r   = opt_res_.z(acc_num + 1) / x0;
	  ofs << "opt_et_basis: "
	    << n   << " " 
	    << num << " "
	    << x0 << " "
	    << r << endl;
	  acc_num += num;
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

      PrintInDebug("keys_values read data");
      keys_values_.Read(ifs);
      
      PrintInDebug("check data type");
      try {
	this->convertData();
      } catch (exception& e) {
	string msg; SUB_LOCATION(msg);
	msg += "failed to convertData. ";
	msg += "Error message is: \n";
	msg += e.what();
	throw runtime_error(msg);
      }

      PrintInDebug("create factory");
      IFactory* ptr = CreateFactory(keys_values_);
      factory_ = boost::shared_ptr<IFactory>(ptr);
      
      PrintInDebug("setting opt target");
      try {
	this->setOptTarget();
      } catch(exception& e) {
	string msg; SUB_LOCATION(msg);
	msg += "\nfailed when setting opt targert\n";
	msg += e.what();
	throw runtime_error(msg);
      }

      PrintInDebug("setting optimizer");
      this->setOptimizer();
          
      timer_.End("read");
    }
    void Compute() {

      timer_.Start("calc");

      if(optimizer_.get()){}else
	err_OptimizerNull();

      if(opt_target_.get()){}else
	err_OptTargetNull();

      if(zs_.rows() == 0) 
	err_zsNull();

      PrintInDebug("Prints objects");
      if(debug_lvl_ > 0) {
	optimizer_->Print();
	opt_target_->Display();
	cout << "zs.rows: " << zs_.rows() << endl;
      }

      try {
	IOptimizer<CD>::Func func = 
	  bind(&IOptTarget::Compute, opt_target_,
	       _1, _2, _3, _4);
	opt_res_ = optimizer_->Optimize(func, zs_);

      } catch (exception& e) {
	cout << "Error on computing ";
	cout << e.what() << endl;
      }
      
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
  OptCBFController::OptCBFController() : impl_(new Impl(0)) {}
  OptCBFController::OptCBFController(int d) : impl_(new Impl(d)) {}
  OptCBFController::~OptCBFController() { 
    delete impl_; }
  void OptCBFController::Read(const char* fn) { 
    impl_->Read(fn); }
  void OptCBFController::Compute() { impl_->Compute(); }
  void OptCBFController::Write(const char* fn) { 
    impl_->Write(fn); }
}

