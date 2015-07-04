
#include <stdexcept>
#include <iostream>
#include "controller.hpp"

using namespace std;
using namespace opt_cbf_h;

void PrintExample() {

  cout << "channel: 1s->kp\n";
  cout << "dipole: length\n";
  cout << "energy: 0.6\n";
  cout << "basis_type: STO\n";
  cout << "opt_basis: 2 (0.8,-0.1)\n";
  cout << "opt_basis: 2 (0.4,-0.6)\n";
  cout << "write_hess: true hess.out\n";
  cout << "write_psi:  true psi.out\n";
  cout << "max_iter: 50\n";
  cout << "eps: 0.00001\n";
  
}

void PrintHelp() {
  cout << "usage: opt_cbf input_file output_file [calculation option]" << endl;
  cout << "       opt_cbf (print option)" << endl;
  cout << "read input_file and compute. " << endl;
  cout << "optimization of orbital exponents of CBF." << endl;
  cout << endl;
  cout << "calculation option are" << endl;
  cout << "-d, --debug           : switch to debugging mode" << endl;
  cout << endl;
  cout << "print option are" << endl;
  cout << "-h, --help            : print help" << endl;
  cout << "-e, --example         : print example input file" << endl;
  cout << endl;
}

int main(int argc, char* argv[]) {

  if(argc == 1) {
    PrintHelp();
    return 0;
  }

  if(argc == 2) {
    if(string(argv[1]) == "-h" || string(argv[1]) == "--help") 
      PrintHelp();
    else if(string(argv[1]) == "-e" || string(argv[1]) == "--example")
      PrintExample();
    else
      PrintHelp();
    return 0;
  }

  int debug_lvl(0);
  if(argc >= 4) {
    string option_str(argv[3]);
    if(option_str == "-d")
      debug_lvl = 1;
  }

  OptCBFController controller(debug_lvl);
  char* in_file = argv[1];
  char* out_file= argv[2];

  if(debug_lvl > 0)
    cout << "Reading file" << endl;
  
  try {
    controller.Read(in_file);
  } catch (exception& e) {
    cout << "Error on reading file" << endl;
    cout << "Error message :" << endl;
    cout << e.what() << endl;
    return 1;
  }

  if(debug_lvl > 0)
    cout << "Computing" << endl;  

  try {
    controller.Compute();
  } catch (exception& e) {
    cout << "Error on computing " << endl;
    cout << "Error message : " << endl;
    cout << e.what() << endl;
    return 1;
  }

  if(debug_lvl > 0)
    cout << "Writing" << endl;  

  try {
    controller.Write(out_file);
  } catch (exception& e) {
    cout << "Error on writing file." << endl;
    cout << "Error message : " << endl;
    cout << e.what() << endl;
    return 1;
  }
}


