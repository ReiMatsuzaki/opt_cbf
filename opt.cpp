#include <iostream>
#include <Eigen/Core>
#include <Eigen/LU>
#include "opt.hpp"

namespace {
  using namespace Eigen;
  using namespace std;
}

namespace opt_cbf_h {

  // ============= Interface ========================
  template<class F>
  IOptimizer<F>::~IOptimizer() {}
  
  // ============= Newton with Restriction ==========
  template<class F>
  OptimizerRestricted<F>::OptimizerRestricted
  (int _max_iter, double _eps, IRestriction<F>* ptr) :
    max_iter_(_max_iter), eps_(_eps), 
    restriction_(ptr), debug_level_(0) {}

  template<class F>
  OptimizerRestricted<F>::OptimizerRestricted
  (int _max_iter, double _eps, IRestriction<F>* ptr, int _d) :
    max_iter_(_max_iter), eps_(_eps), 
    restriction_(ptr), debug_level_(_d) {}

  template<class F>
  OptimizerRestricted<F>::~OptimizerRestricted() {
    delete restriction_;
  }

  template<class F>
  void PrintDebug(OptRes<F> opt_res) {
    cout << opt_res.z << endl;
    cout << opt_res.grad << endl;
  }
  template<class F>
  OptRes<F> OptimizerRestricted<F>::Optimize
  (FuncValGradHess f, VecF z0){

    int num = z0.rows();

    // initialize
    OptRes<F> res;
    res.convergence = false;
    res.z = z0;
    res.value= F(0);
    res.hess = MatF::Zero(num, num);
    res.grad = VecF::Zero(num);
    VecF dz_full   = VecF::Zero(num);

    int n1 = restriction_->size();
    VecF grad(n1);
    MatF hess(n1, n1);
    VecF dz(n1);

    // start loop
    for(int i = 0; i < max_iter_; i++) {

      // compute grad and hess
      f(res.z, &res.value, &res.grad, &res.hess);

      // restriction
      restriction_->SetVars(res.z);
      grad = restriction_->Grad(res.grad);
      hess = restriction_->Hess(res.grad, res.hess);

      // print if debugging mode
      if(debug_level_ > 0)
	PrintDebug(res);
      
      // update vector in restriction space
      dz = -hess.fullPivLu().solve(grad);

      // update in normal space
      restriction_->Shift(dz);
      res.z = restriction_->Xs();
      res.iter_num = i + 1;

      // check convergence
      double dz_norm = std::abs(dz.norm());
      bool check1 = dz_norm < eps_;
      double max_grad = res.grad.array().abs().maxCoeff();
      bool check2 =  max_grad < eps_;
      if( check1 && check2) {
	res.convergence = true;
	break;
      }
    }
    
    return res;
    
  }
  
  // ============= Simple Newton ====================
  template<class F>
  OptimizerNewton<F>::OptimizerNewton(int _max_iter, double _eps) {
    IRestriction<F>* no_rist = new NoRestriction<F>();
    optimizer_ = new OptimizerRestricted<F>(_max_iter, _eps, no_rist);
  }
  template<class F>
  OptimizerNewton<F>::OptimizerNewton(int _max_iter, double _eps,
				      int _d_lvl) {
    IRestriction<F>* no_rist = new NoRestriction<F>();
    optimizer_ = new OptimizerRestricted<F>(_max_iter, _eps, no_rist, _d_lvl);    
  }
  template<class F>
  OptimizerNewton<F>::~OptimizerNewton() {
    delete optimizer_;
  }
  template<class F>
  OptRes<F> OptimizerNewton<F>::Optimize(FuncValGradHess f,VecF z0) {
    return optimizer_->Optimize(f, z0);
  }

  // ============== explicit instance ===============
  typedef std::complex<double> CD;
  template class OptimizerNewton<double>;
  template class OptimizerNewton<CD>;
  template class OptimizerRestricted<double>;
  template class OptimizerRestricted<CD>;
  
}



