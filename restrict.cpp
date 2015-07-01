#include <iostream>
#include <stdexcept>
#include <boost/lexical_cast.hpp>
#include <Eigen/Core>
#include <macros.hpp>

#include "restrict.hpp"

namespace {
  using std::cout;
  using std::endl;
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
    VecF xs = VecF::Zero(this->size());

    int ni_m1(0);
    for(CIT it = num_x0_r_list_.begin(),
	  it_end = num_x0_r_list_.end();
	it != it_end; ++it) {
      
      int ni = ni_m1 + get<0>(*it);
      F   a   = get<1>(*it);
      F   r   = get<2>(*it);

      for(int k = ni_m1; k < ni; k++) {
	xs(k) = a * pow(r, k - ni_m1);
      }

      ni_m1 = ni;
    }

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

    // see ipython notebook for variable notations
    int ni_m1(0); // represent n_{i-1}
    for(int i = 0; i < num_rest; i++) {

      F df_da(0), df_dr(0);
      int ki_max = get<0>(num_x0_r_list_[i]);
      F   ai     = get<1>(num_x0_r_list_[i]);
      F   ri     = get<2>(num_x0_r_list_[i]);
      int ni = ni_m1 + ki_max;
      for(int k = ni_m1; k < ni; k++) {
	df_da += g(k) * pow(ri, k - ni_m1);
	df_dr += g(k) * ai * F(k-ni_m1) * pow(ri, k-ni_m1-1);
      }
      grad(2*i)   = df_da;
      grad(2*i+1) = df_dr;
      ni_m1 = ni;
    }
    
    return grad;
  }
  template<class F>
  typename MultiEvenTemp<F>::MatF  MultiEvenTemp<F>::Hess
  (const VecF& g, const MatF& h) const {

    if(g.rows() != this->size() ||
       g.rows() != h.rows() ||
       g.rows() != h.cols() ) {
      std::string msg; SUB_LOCATION(msg);
      msg += "\nsize is invalid\n"; 
      throw std::runtime_error(msg);
    }

    int num_rest = num_x0_r_list_.size();
    MatF hess(num_rest * 2, num_rest * 2);

    int ni_m1(0);
    for(int i = 0; i < num_rest; i++) {

      int ki_max = get<0>(num_x0_r_list_[i]);
      F   ai     = get<1>(num_x0_r_list_[i]);
      F   ri     = get<2>(num_x0_r_list_[i]);
      int ni = ni_m1 + ki_max;

      int nj_m1(0);
      for(int j = 0; j < num_rest; j++) {

	int kj_max = get<0>(num_x0_r_list_[j]);
	F   aj     = get<1>(num_x0_r_list_[j]);
	F   rj     = get<2>(num_x0_r_list_[j]);
	int nj     = nj_m1 + kj_max;
	
	F d2f_dada(0), d2f_dadr(0), d2f_drdr(0);
	for(int k = ni_m1; k < ni; k++) {
	  if(i == j) {
	    F gk = g(k);
	    d2f_dadr += gk * F(k-ni_m1) * pow(ri, k-ni_m1-1);
	    d2f_drdr += gk * ai * F((k-ni_m1)*(k-ni_m1-1)) * pow(ri, k-ni_m1-2);
	  }
	  for(int l = nj_m1; l < nj; l++) {
	    F hkl = h(k, l);
	    d2f_dada += hkl * pow(ri, k - ni_m1) * pow(rj, l - nj_m1);
	    d2f_dadr += hkl * pow(ri, k - ni_m1) * 
	      aj * F(l-nj_m1) * pow(rj, l-nj_m1-1);
	    d2f_drdr += hkl * 
	      ai * F(k-ni_m1) * pow(ri, k-ni_m1-1) * 
	      aj * F(l-nj_m1) * pow(rj, l-nj_m1-1);
	  }
	}
	nj_m1 = nj;

	hess(2*i,   2*j  ) = d2f_dada;
	hess(2*i+1, 2*j  ) = d2f_dadr;
	hess(2*i,   2*j+1) = d2f_dadr;
	hess(2*i+1, 2*j+1) = d2f_drdr;
      }
      ni_m1 = ni;
    }

    return hess;
  }
  template<class F>
  void MultiEvenTemp<F>::Shift(const VecF& dz)  { 
    
    if(num_x0_r_list_.size() * 2 != dz.rows()) {
      string msg; SUB_LOCATION(msg);
      msg += "\ninvalid size\n";
      throw std::runtime_error(msg);
    }

    int i(0);
    for(IT it = num_x0_r_list_.begin(),
	  it_end = num_x0_r_list_.end();
	it != it_end; ++it) {
      int m = get<0>(*it);
      F   a = get<1>(*it);
      F   r = get<2>(*it);
      *it = make_tuple(m, a + dz(2*i), r + dz(2*i+1));
      i++;
    }
  }

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
