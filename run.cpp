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
  cout << "usage: opt_cbf input_file" << endl;
  cout << "       opt_cbf option" << endl;
  cout << "read input_file and compute. " << endl;
  cout << "optimization of orbital exponents of CBF." << endl;
  cout << endl;
  cout << "-h, --help            : print help" << endl;
  cout << "-e, --example         : print example input file";
  cout << endl;
}

int main(int argc, char* argv[]) {

  if(argc == 2) {
    if(string(argv[1]) == "-h" || string(argv[1]) == "--help") 
      PrintHelp();
    else if(string(argv[1]) == "-e" || string(argv[1]) == "--example")
      PrintExample();
    return 0;
  }

  if(argc != 3) {
    PrintHelp();
    return 0;
  }

  OptCBFController controller;
  char* in_file = argv[1];
  char* out_file= argv[2];  
  try {
    controller.Read(in_file);
  } catch (exception& e) {
    cout << "Error on reading file" << endl;
    cout << "Error message :" << endl;
    cout << e.what() << endl;
    exit(1);
  }

  try {
    controller.Compute();
  } catch (exception& e) {
    cout << "Error on computing " << endl;
    cout << "Error message : " << endl;
    cout << e.what() << endl;
    exit(1);
  }

  try {
    controller.Write(out_file);
  } catch (exception& e) {
    cout << "Error on writing file." << endl;
    cout << "Error message : " << endl;
    cout << e.what() << endl;
    exit(1);
  }
}


