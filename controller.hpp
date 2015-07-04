#ifndef CONTROLLER_HPP
#define CONTROLLER_HPP

#include <string>

namespace {
  using std::string;
}

namespace opt_cbf_h {

  // optimization of CBF
  // this class is designed using IMPL idiom
  class OptCBFController {
  private:
    // -------- Member Field ---------
    class Impl;
    Impl* impl_;

    // -------- Stop copy ------------
    OptCBFController(OptCBFController&);
    
  public:
    // -------- Constructor ----------
    OptCBFController();
    explicit OptCBFController(int _debug_lvl);
    ~OptCBFController();
    // -------- Method ---------------
    void Read(const char* filename);
    void Compute();
    void Write(const char* filename);
  };
}

#endif
