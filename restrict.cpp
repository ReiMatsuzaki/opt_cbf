#include <iostream>
#include <stdexcept>
#include <boost/lexical_cast.hpp>
#include <Eigen/Core>
#include <macros.hpp>

#include "restrict.hpp"

namespace {
  using boost::get;
  using std::string;
  using boost::lexical_cast;
}

namespace opt_cbf_h {

  // ========== Interface =================

  template<class F>
  IRestriction<F>::~IRestriction() {}

  // ========== EvenTemp =================
  template<class F>
  EvenTemp<F>::EvenTemp() : num_(1), x0_(1.0), ratio_(1.0) {}
  template<class F>
  EvenTemp<F>::~EvenTemp() {}
  template<class F>
  void EvenTemp<F>::SetVars(const Matrix<F, Dynamic, 1>& xs) {

    num_   = xs.rows();
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
  Matrix<F, Dynamic, 1> EvenTemp<F>::Grad(const Matrix<F, Dynamic, 1>& dxs) const {

    int num = dxs.rows();
    
    F df_da = F(0);
    F df_dr = F(0);
    for(int i = 0; i < num; i++) {
      df_da += pow(ratio_, i) * dxs(i);
      df_dr += x0_ * F(i) * pow(ratio_, i-1) * dxs(i);
    }

    VecF grad(2);
    grad(0) = df_da; grad(1) = df_dr;

    return grad;

  }
  /**
   * represent even-tempered sequence
   * x_n = a r^n (n=0,...,N-1)
   * df/dada   = sum_nm dx_n/da dx_m/da df^2/(dx_n)(dx_m) 
   *           = sum_nm r^n r^m f_nm
   * d^2f/dadr = d_a sum_n anr^{n-1} f_n
   *           = sum_n nr^{n-1}f_n + sum_nm anr^{n+m-1} f_nm
   * d^2f/drdr = d_r sum_n anr^n{n-1} f_n
   *           = an(n-1)r^{n-2} f_n
   *             + n a r^{n-1} m a r^{m-1} f_nm
   */
  template<class F>
  Matrix<F, Dynamic, Dynamic> EvenTemp<F>::Hess
  (const Matrix<F, Dynamic, 1>&       dxs, 
   const Matrix<F, Dynamic, Dynamic>& ddxs) const {
						
    int num = ddxs.rows();
    F d2f_dada = F(0);
    F d2f_dadr = F(0);
    F d2f_drdr = F(0);

    for(int n = 0; n < num; n++) {
      d2f_dadr += F(n) * pow(ratio_, n-1) * dxs(n);
      d2f_drdr += F(n*(n-1)) * x0_ * pow(ratio_, n-2) * dxs(n);
      for(int m = 0; m < num; m++) {
	F h = ddxs(n, m);
	d2f_dada += pow(ratio_, n+m)   * h;
	d2f_dadr += pow(ratio_, n+m-1) * F(m) * x0_ * h;
	if(m != 0 && n != 0)
	  d2f_drdr += pow(ratio_, n+m-2) * F(n*m) *x0_*x0_ * h;
      }
    }
      
    Matrix<F, Dynamic, Dynamic> hess(2, 2);
    hess << 
      d2f_dada, d2f_dadr, 
      d2f_dadr, d2f_drdr;
    return hess;
  }
  template<class F>
  Matrix<F, Dynamic, 1> EvenTemp<F>::Xs() const {

    VecF xs(num_);
    F val = x0_;
    for( int n = 0; n < num_; n++) {
      xs(n) = val;
      val *= ratio_;
    }

    return xs;
  }
  template<class F>
  void EvenTemp<F>::Shift(const VecF& ys) {
    x0_ += ys(0);
    ratio_ += ys(1);
  }

  // ============ Multi Even Temp ====================
  template<class F>
  MultiEvenTemp<F>::MultiEvenTemp
  (const vector<int>& is) {
   
    int num = is.size();
    num_x0_r_list_.resize(num);

    for(int n = 0; n < num; n++) 
      num_x0_r_list_[n] = make_tuple(is[n], F(0), F(0));

  }
  template<class F>
  MultiEvenTemp<F>::~MultiEvenTemp() {}
  template<class F>
  typename MultiEvenTemp<F>::IFF
  MultiEvenTemp<F>::num_x0_r(int i) const {
    
    if(i >= num_x0_r_list_.size()) {
      std::string msg; SUB_LOCATION(msg);
      msg += "\n index exceed size of vector\n";
      throw std::runtime_error(msg);
    }

    return num_x0_r_list_[i];
  }
  template<class F>
  void MultiEvenTemp<F>::SetVars(const VecF& xs) {

    int acc(0);
    for(IT it = num_x0_r_list_.begin(),
	  it_end = num_x0_r_list_.end();
	it != it_end; ++it) {

      if(acc+1 >= xs.size() ) {
	std::string msg; SUB_LOCATION(msg);
	msg += "\n acc exceed xs.size \n";
	throw std::runtime_error(msg);
      }

      int num = get<0>(*it);
      F   x0  = xs(acc);
      F   r   = xs(acc + 1) / xs(acc);
      acc += num;
      *it = make_tuple(num, x0, r);
    }
  }
  template<class F> typename MultiEvenTemp<F>::VecF 
  MultiEvenTemp<F>::Xs() const {
    VecF xs(1);
    return xs;
  }

  template<class F>
  int MultiEvenTemp<F>::size() const { 

    int acc(0);
    for(CIT it = num_x0_r_list_.begin(),
	  it_end = num_x0_r_list_.end(); 
	it != it_end; ++it) {
      acc += get<0>(*it);
    }
    return acc;
  }
  template<class F>
  typename MultiEvenTemp<F>::VecF
  MultiEvenTemp<F>::Grad(const VecF& g) const {

    if(g.rows() != this->size()) {
      std::string msg; SUB_LOCATION(msg);
      msg += "\nsize is invalid\n";
      msg += boost::lexical_cast<string>(g.rows());
      msg += " ";
      msg += boost::lexical_cast<string>(this->size());
    }

    int num_rest = num_x0_r_list_.size();
    VecF grad(num_rest * 2);

    int n0(0);
    for(int i = 0; i < num_rest; i++) {

      F df_da(0);
      F df_dr(0);
      int num = get<0>(num_x0_r_list_[i]);
      F   x0  = get<1>(num_x0_r_list_[i]);
      F    r  = get<2>(num_x0_r_list_[i]);
      for(int n = 0; n < num; n++) {
	df_da += pow(r, n) * g(n0 + n);
	df_dr += x0 * F(n) * pow(r, n-1) * g(n0 + n);
      }
      grad(2*i)   = df_da;
      grad(2*i+1) = df_dr;
	
      n0 += num;
    }

    return grad;
  }
  template<class F>
  typename MultiEvenTemp<F>::MatF  MultiEvenTemp<F>::Hess
  (const VecF& g, const MatF& h) const {
    return MatF(h);
  }
  template<class F>
  void MultiEvenTemp<F>::Shift(const VecF& dz)  { }

  // ============ Explicit Instance ==================
  typedef std::complex<double> CD;
  template class IRestriction<double>;
  template class IRestriction<CD>;
  template class NoRestriction<double>;
  template class NoRestriction<CD>;
  template class EvenTemp<double>;
  template class EvenTemp<CD>;
  template class MultiEvenTemp<double>;
  template class MultiEvenTemp<CD>;
  
}
