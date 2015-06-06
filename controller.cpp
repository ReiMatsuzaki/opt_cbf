#include <iostream>
#include <boost/scoped_ptr.hpp>
#include <keys_values.hpp>
#include "controller.hpp"
#include "opt_cbf.hpp"
#include "opt.hpp"

namespace {
  using std::cout;
  using std::endl;
  using std::ifstream;
  typedef std::complex<double> CD;
  using boost::scoped_ptr;
}

namespace opt_cbf_h {

  class OptCBFController::Impl {
  public:
    // -------- Member Field ---------
    KeysValues keys_values_;
    scoped_ptr<IOptTarget> opt_target_;
    scoped_ptr<IOptimizer<CD> > optimizer_;

    // -------- Constructor ---------
    Impl() : keys_values_(":", " ") {}

    // ------- Read support ---------
    void convertData() {

      keys_values_.ConvertValues<string>("basis");
      keys_values_.ConvertValues<int>("max_iter");
      keys_values_.ConvertValues<double>("eps");
      
    }
    void setOptTarget() {

    }
    void setOptimizer() {
      
    }
    
    // -------- Method --------------
    void Read(const char* filename) {
      
      ifstream ifs(filename);

      if(ifs.fail()) {
	cout << "Failed to open file." << endl;
	cout << "File name: " << filename << endl;
	exit(1);
      }
      
      keys_values_.Read(ifs);

      this->convertData();
      this->setOptTarget();
      this->setOptimizer();
      
    }
    void Compute() {

    }
    void Write(const char* filename) {

    }


    
  };

  // =========== interface =============
  OptCBFController::OptCBFController() : impl_(new Impl()) {}
  OptCBFController::~OptCBFController() { delete impl_; }
  void OptCBFController::Read(const char* fn) { impl_->Read(fn); }
  void OptCBFController::Compute() { impl_->Compute(); }
  void OptCBFController::Write(const char* fn) { impl_->Write(fn); }
}
