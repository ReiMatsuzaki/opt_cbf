#ifndef OPT_HPP
#define OPT_HPP

#include <Eigen/Core>
#include <tr1/functional>

namespace {
  using namespace Eigen;
  using std::tr1::function;
}

namespace opt_cbf_h {

  // result of optimization
  template<class F>
  class OptRes {
  public:
    OptRes() :
      convergence(false), iter_num(0){}
    bool convergence;
    Matrix<F, Dynamic, 1> z;
    int iter_num;
  };

  // interface for optimization class
  template<class F>
  class IOptimizer {
    typedef Matrix<F, Dynamic, 1>       VecF;
    typedef Matrix<F, Dynamic, Dynamic> MatF;
    //    typedef function<void (const VecF&, VecF*, MatF*)> FuncGradHess;
    typedef function<void (const VecF&, F*, VecF*, MatF*)> FuncValGradHess;
  public:
    virtual ~IOptimizer() {}
    //    virtual OptRes<F> Optimize(FuncGradHess f, VecF z0) = 0;
    virtual OptRes<F> Optimize(FuncValGradHess f, VecF z0) = 0;
  };

  // simple newton iteration
  template<class F>
  class OptimizerNewton : public IOptimizer<F> {
    typedef Matrix<F, Dynamic, 1>       VecF;
    typedef Matrix<F, Dynamic, Dynamic> MatF;
    //    typedef function<void (VecF, VecF*, MatF*)> FuncGradHess;
    typedef function<void (const VecF&, F*, VecF*, MatF*)> FuncValGradHess;
  public:
    OptRes<F> Optimize(FuncValGradHess f, VecF z0) {

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
	dz = hess.colPivHouseholderQr().solve(grad);
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
  };
}

#endif
