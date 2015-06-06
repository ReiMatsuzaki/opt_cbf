#include <Eigen/Core>
#include <Eigen/LU>
#include "opt.hpp"

namespace {
  using namespace Eigen;
}

namespace opt_cbf_h {

  template<class F>
  OptRes<F> OptimizerNewton<F>::Optimize
  (FuncValGradHess f,VecF z0) {

      // initialize
      OptRes<F> res;
      res.convergence = false;

      res.z = VecF::Zero(2);
      res.z = z0;
      double eps = 0.000000001;
      int max_num = 100;

      int num = z0.rows();
      F    value;
      MatF hess = MatF::Zero(num, num);
      VecF grad = VecF::Zero(num);
      VecF dz = VecF::Zero(num);

      // start loop
      for(int i = 0; i < max_num; i++) {

	// compute grad and hess
	f(res.z, &value, &grad, &hess);

	// update
	dz = hess.fullPivLu().solve(grad);
	res.z -= dz;

	// check convergence
	double dz_norm = std::abs(dz.norm());
	bool check1 = dz_norm < eps;
	double max_grad = grad.array().abs().maxCoeff();
	bool check2 =  max_grad < eps;
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



