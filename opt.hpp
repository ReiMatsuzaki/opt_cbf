#ifndef OPT_HPP
#define OPT_HPP

#include <Eigen/Core>
#include <boost/function.hpp>
//#include <boost/any.hpp>

namespace {
  using std::string;
  //  using std::map;
  using namespace Eigen;
  using boost::function;
  //  using boost::any;
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
    F                     value;
    Matrix<F, Dynamic, 1> grad;
    Matrix<F, Dynamic, Dynamic> hess;

    int iter_num;
    //    map<string, any> others_map;
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
  private:
    // --------------- type ----------------
    typedef Matrix<F, Dynamic, 1>       VecF;
    typedef Matrix<F, Dynamic, Dynamic> MatF;
    typedef function<void (const VecF&, F*, VecF*, MatF*)>
    FuncValGradHess;

    // ------------- Field -----------------
    int max_iter_;
    double eps_;
    int debug_level_;

  public:
    OptimizerNewton(int _max_iter, double _eps);
    OptimizerNewton(int _max_iter, double _eps, int _d_lvl);
    OptRes<F> Optimize(FuncValGradHess f, VecF z0);
  };

  /*
  // Decorator pattern
  template<class F>
  class Decorator : public IOptimizer<F> {
  private:
    // --------------- type ----------------
    typedef Matrix<F, Dynamic, 1>       VecF;
    typedef Matrix<F, Dynamic, Dynamic> MatF;
    typedef function<void (const VecF&, F*, VecF*, MatF*)>
    FuncValGradHess;

    // --------------- Field ---------------
    IOptimizer* optimizer_;
  };
  */
}

#endif
