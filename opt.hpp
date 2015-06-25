#ifndef OPT_HPP
#define OPT_HPP

#include <Eigen/Core>
#include <boost/function.hpp>
#include "restrict.hpp"

namespace {

  using std::string;
  using namespace Eigen;
  using boost::function;

}

namespace opt_cbf_h {

  //-------- forward declaration ---------------
  // template<class F> class IRestriction;
  
  // ========== result of optimization =============
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

  // ========== interface for optimization class =====
  template<class F>
  class IOptimizer {
    typedef Matrix<F, Dynamic, 1>       VecF;
    typedef Matrix<F, Dynamic, Dynamic> MatF;
    typedef function<void (const VecF&, F*, VecF*, MatF*)>
    FuncValGradHess;
  public:
    virtual ~IOptimizer();
    virtual OptRes<F> Optimize(FuncValGradHess f, VecF z0) = 0;
  };

  // =========== with restrictions ==================
  template<class F>
  class OptimizerRestricted : public IOptimizer<F> {
  private:
    // ------------ Type --------------------
    typedef Matrix<F, Dynamic, 1>       VecF;
    typedef Matrix<F, Dynamic, Dynamic> MatF;
    typedef function<void (const VecF&, F*, VecF*, MatF*)>
    FuncValGradHess;
    // ------------ Field -------------------
    int max_iter_;
    double eps_;
    IRestriction<F>* restriction_;
    int debug_level_;
    // ------------ Uncopyable ---------------
    OptimizerRestricted(const OptimizerRestricted<F>&);
    
  public:
    // ------------ Constructors ------------
    OptimizerRestricted(int, double, IRestriction<F>*);
    OptimizerRestricted(int, double, IRestriction<F>*, int);
    ~OptimizerRestricted();
    
    // ------------ Method ------------------
    OptRes<F> Optimize(FuncValGradHess f, VecF z0);
  };

  // =========== simple newton iteration ============
  template<class F>
  class OptimizerNewton : public IOptimizer<F> {
  private:
    // --------------- type ----------------
    typedef Matrix<F, Dynamic, 1>       VecF;
    typedef Matrix<F, Dynamic, Dynamic> MatF;
    typedef function<void (const VecF&, F*, VecF*, MatF*)>
    FuncValGradHess;

    // ------------- Field -----------------
    OptimizerRestricted<F>* optimizer_;

    // ------------ Uncopyable -------------
    OptimizerNewton(const OptimizerNewton&);

  public:
    // ------------ Constructors ------------
    OptimizerNewton(int _max_iter, double _eps);
    OptimizerNewton(int _max_iter, double _eps, int _d_lvl);
    ~OptimizerNewton();
    OptRes<F> Optimize(FuncValGradHess f, VecF z0);
  };
}

#endif
