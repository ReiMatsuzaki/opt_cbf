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
    bool convergence;
    Matrix<F, Dynamic, 1> z;
  };

  // interface for optimization class
  template<class F>
  class IOptimizer {
    typedef Matrix<F, Dynamic, 1>       VecF;
    typedef Matrix<F, Dynamic, Dynamic> MatF;
    typedef function<void (VecF, VecF*, MatF*)> FuncGradHess;
  public:
    virtual ~IOptimizer() {}
    virtual OptRes<F> Optimize(FuncGradHess f, VecF z0) = 0;
  };

  // simple newton iteration
  template<class F>
  class OptByNewton : public IOptimizer<F> {
    typedef Matrix<F, Dynamic, 1>       VecF;
    typedef Matrix<F, Dynamic, Dynamic> MatF;
    typedef function<void (VecF, VecF*, MatF*)> FuncGradHess;
  public:
    OptRes<F> Optimize(FuncGradHess f, VecF z0) {

      // initialize
      OptRes<F> res;
      res.convergence = false;

      res.z = VecF::Zero(2);
      res.z = z0;
      double eps = 0.0000001;
      int max_num = 20;


      int num = z0.rows();
      MatF hess = MatF::Zero(num, num);
      VecF grad = VecF::Zero(num);
      VecF dz = VecF::Zero(num);

      // start loop
      for(int i = 0; i < max_num; i++) {

	// compute grad and hess
	f(res.z, &grad, &hess);

	// update
	dz = hess.colPivHouseholderQr().solve(grad);
	res.z -= dz;

	// check convergence
	if( dz.norm() < eps) {
	  res.convergence = true;
	  break;
	}
      }
      return res;
    }
  };

}

#endif
