#include "restrict.hpp"
#include <Eigen/Core>

namespace opt_cbf {

  template<class F>
  void EvenTemp<F>::SetVars(const VecF& xs) {

    x0_    = xs(0);
    ratio_ = xs(1) / xs(0);

  }

  /**
   * represent eventempered sequence
   * x_n = a r^n (n = 0, ..., N-1)
   * df/da = sum_n (dx_n/da)(df/dx_n)
   *       = sum_n r^n df/dx_n
   * df/dr = sum_n (dx_n/dr)(df/dx_n)
   *       = sum_n anr^{n-1} df/dx_n
   */      
  template<class F>
  VecF EvenTemp<F>::Grad(const VecF& dxs) {

    int num = dxs.rows();
    
    F df_da = F(0);
    F df_dr = F(0);
    for(int i = 0; i < num; i++) {
      df_da += pow(ratio_, i) * dxs(0);
      df_dr += x0_ * F(i) * pow(ratio_, i) * dxs(0);
    }

    VecF grad(2);
    grad(0) = df_da; grad(1) = df_dr;

    return grad(0);
  }

}
