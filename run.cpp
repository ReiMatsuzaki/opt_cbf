#include "controller.hpp"

using namespace opt_cbf_h;

int main() {

  OptCBFController controller;
  controller.Read("supply/sample.in");
  controller.Compute();
  controller.Write("supply/sample.out");
  
}


