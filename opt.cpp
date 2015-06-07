#include <iostream>
#include <Eigen/Core>
#include <Eigen/LU>
#include "opt.hpp"

namespace {
  using namespace Eigen;
  using namespace std;
}

namespace opt_cbf_h {
  
  template<class F>
  OptimizerNewton<F>::OptimizerNewton
  (int _max_iter, double _eps) :
    max_iter_(_max_iter), eps_(_eps), debug_level_(0) {}
  template<class F>
  OptimizerNewton<F>::OptimizerNewton
  (int _max_iter, double _eps, int _d_lvl) :
    max_iter_(_max_iter), eps_(_eps), debug_level_(_d_lvl) {}

  template<class F>
  OptRes<F> OptimizerNewton<F>::Optimize
  (FuncValGradHess f,VecF z0) {

    int num = z0.rows();

    // initialize
    OptRes<F> res;
    res.convergence = false;

    res.z = z0;
      
    F    value;
    MatF hess = MatF::Zero(num, num);
    VecF grad = VecF::Zero(num);
    VecF dz   = VecF::Zero(num);

    // start loop
    for(int i = 0; i < max_iter_; i++) {
      
      // compute grad and hess
      f(res.z, &value, &grad, &hess);

      // print if debugging mode
      if (debug_level_ > 0) {
	cout << i << "; ";
	for(int i = 0; i < num; i++)
	  cout << res.z(i) << ", ";
	cout << "; ";
	for(int i = 0; i < num; i++)
	  cout << grad(i) << ", ";
	cout << endl;
      }
      
      // update
      dz = hess.fullPivLu().solve(grad);
      res.z -= dz;
      res.iter_num = i + 1;

      // check convergence
      double dz_norm = std::abs(dz.norm());
      bool check1 = dz_norm < eps_;
      double max_grad = grad.array().abs().maxCoeff();
      bool check2 =  max_grad < eps_;
      if( check1 && check2) {
	res.convergence = true;
	break;
      }
    }
    return res;
  }

  // explicit instance
  template class OptimizerNewton<double>;
  template class OptimizerNewton<std::complex<double> >;
}



