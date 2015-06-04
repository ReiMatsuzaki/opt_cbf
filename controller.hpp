#ifndef CONTROLLER_HPP
#define CONTROLLER_HPP

#include "opt_cbf.hpp"
#include "opt.hpp"

namespace {
  using std::string;
}

namespace opt_cbf_h {
  class OptCBFController {
  private:
    // interface to access other basis
    
  public:
    // method
    void Read(const string& filename);
    void Compute();
    void Write(const string& filename);
  };
}

#endif
