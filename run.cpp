#include <stdexcept>
#include <iostream>
#include "controller.hpp"

using namespace std;
using namespace opt_cbf_h;

int main() {

  OptCBFController controller;
  try {
    controller.Read("supply/sample.in");
  } catch (exception& e) {
    cout << "Error on reading file" << endl;
    cout << "Error message :" << endl;
    cout << e.what() << endl;
  }

  try {
    controller.Compute();
  } catch (exception& e) {
    cout << "Error on computing " << endl;
    cout << "Error message : " << endl;
    cout << e.what() << endl;
  }

  try {
    controller.Write("supply/sample.out");
  } catch (exception& e) {
    cout << "Error on writing file." << endl;
    cout << "Error message : " << endl;
    cout << e.what() << endl;
  }
}


