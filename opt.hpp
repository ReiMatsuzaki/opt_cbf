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
    typedef function<void (const VecF&, F*, VecF*, MatF*)>
    FuncValGradHess;
  public:
    virtual ~IOptimizer() {}
    virtual OptRes<F> Optimize(FuncValGradHess f, VecF z0) = 0;
  };

  // simple newton iteration
  template<class F>
  class OptimizerNewton : public IOptimizer<F> {
    typedef Matrix<F, Dynamic, 1>       VecF;
    typedef Matrix<F, Dynamic, Dynamic> MatF;
    typedef function<void (const VecF&, F*, VecF*, MatF*)>
    FuncValGradHess;
  public:
    OptRes<F> Optimize(FuncValGradHess f, VecF z0);
  };
}

#endif
