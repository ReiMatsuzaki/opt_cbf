#ifndef OPT_HPP
#define OPT_HPP

#include <Eigen/Core>
#include <boost/function.hpp>

namespace {

  using std::string;
  using namespace Eigen;
  using boost::function;

}

namespace opt_cbf_h {

  //-------- forward declaration ---------------
  template<class F> class IRestriction;
  
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

  };

  // ========== interface for optimization class =====
  template<class F>
  class IOptimizer {
    typedef Matrix<F, Dynamic, 1>       VecF;
    typedef Matrix<F, Dynamic, Dynamic> MatF;
    typedef boost::function<void (const VecF&, F*, VecF*, MatF*)>    FuncValGradHess;

  protected:
    // ------------ Field -----------------
    int max_iter_;
    double eps_;
  public:
    // ------------ Constructors ----------
    virtual ~IOptimizer();
    IOptimizer(int m, double e);

    // ------------ Method ----------------
    virtual OptRes<F> Optimize(FuncValGradHess f, VecF z0) = 0;
    // ------------ Accessors --------------
    int max_iter() const { return max_iter_; }
    double eps() const { return eps_; }
  };

  // =========== with restrictions ==================
  template<class F>
  class OptimizerRestricted : public IOptimizer<F> {
  private:
    // ------------ Type --------------------
    typedef Matrix<F, Dynamic, 1>       VecF;
    typedef Matrix<F, Dynamic, Dynamic> MatF;
    typedef boost::function<void (const VecF&, F*, VecF*, MatF*)>
    FuncValGradHess;
    // ------------ Field -------------------
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
    typedef boost::function<void (const VecF&, F*, VecF*, MatF*)>
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
