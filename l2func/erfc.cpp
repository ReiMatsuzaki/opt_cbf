#include "erfc.hpp"

namespace erfc_mori {

  template<class F> bool NormIsLessThan(F x, double y) {
    return std::abs(x) < y;
  }
  double machine_eps() {
    return  2.22 * pow(10.0, -16.0);
  }
  template<class F> int erfc_start_num_term() {
    raise_error<F>();
    return 0; }
  template<> int erfc_start_num_term<double>() {
    return 10; }
  template<> int erfc_start_num_term<CD>() {
    return 10; }

  template<class F>  bool erfc_add_Eh_q(F x, F h) {
    raise_error<F>();
    return true; }
  template<> bool erfc_add_Eh_q(double x, double h) {
    return x < M_PI / h;
  }
  template<> bool erfc_add_Eh_q(CD x, CD h) {
    return x.real() + std::abs(x.imag()) < M_PI / h.real();
  }
  template<class F>
  int ErfcAtN(F x, unsigned int n, F& y, F& yn ) {

    F y0, yn0, one;
    F nnhh, h, xx, t2;

    one = F(1);
    xx = x * x;
    h = sqrt(M_PI / (int)n);

    y0 = F(0);
    for(unsigned int i = 0; i < n; i++) {
      // nnhh = (2.0 * i + 1.0) * h / 2.0;
      nnhh = F(2 * (int)i + 1) * h / 2.0;
      nnhh = nnhh * nnhh;
      t2 = exp(-nnhh) / (nnhh + xx);
      y0 += t2;
      yn0  = t2;
    }

    t2 = F(2) * h / M_PI * exp(-xx) * x;
    y0  *= t2;
    yn0 *= t2;

    if(erfc_add_Eh_q<F>(x, h))
      y0 += 2.0 / (one + exp(2*M_PI*x/h));

    y = y0; yn = yn0;

    return(0);
  }
  template<class F>
  int Erfc(F x, F & y, ErfcCalcData& data) {

    F y0, yn0, ratio;
    double eps = machine_eps();

    unsigned n0, n1;
    n0 = erfc_start_num_term<F>();
    n1 = n0 + 10;

    data.convergence = false;
    for(unsigned n = n0; n < n1; n++) {
      ErfcAtN<F>(x, n, y0, yn0);
      ratio = yn0 / y0;
      data.num_term = n;

      if(data.convergence)
	break;

      if(NormIsLessThan(ratio, eps) &&
	 NormIsLessThan(yn0, eps)) {
	data.convergence = true;
      }
    }

    y = y0;

    if(data.convergence)
      return (0);
    else 
      return (1);
  }
  template<class F>
  int Exp2ErfcAtN(F x, unsigned int n, F& y, F& yn) {

    F y0, yn0;
    F nnhh, h, xx, t2, one;
    xx = x * x;
    h = sqrt(M_PI / (int)n);
    one = F(1);

    y0 = F(0);
    for(unsigned int i = 0; i < n; i++) {
      nnhh = F(2 * (int)i + 1) * h / F(2);
      nnhh = nnhh * nnhh;
      t2 = exp(-nnhh) / (nnhh + xx);
      y0 += t2;
      yn0 = t2;
    }

    t2 = F(2) * h / M_PI * x;
    y0  *= t2;
    yn0 *= t2;

    if(erfc_add_Eh_q<F>(x, h))
      y0 += exp(xx) * F(2) / (one + exp(2*M_PI*x/h));

    y=y0; yn=yn0;
    return(0);
    
  }
  template<class F>
  int Exp2Erfc(F x, F& y, ErfcCalcData& data) {

    F y0, yn0, ratio;
    double eps = machine_eps();

    unsigned n0, n1;
    n0 = erfc_start_num_term<F>();
    n1 = n0 + 10;

    data.convergence = false;
    for(unsigned n = n0; n < n1; n++) {

      Exp2ErfcAtN<F>(x, n, y0, yn0);
      ratio = yn0 / y0;
      data.num_term = n;

      if(data.convergence)
	break;

      if(NormIsLessThan(ratio, eps) &&
	 NormIsLessThan(yn0, eps)) {
	data.convergence = true;
      }
    }

    y = y0;

    if(data.convergence)
      return (0);
    else 
      return (1);    

  }
  
  // explicit instance
  typedef std::complex<double> CD;
  template bool NormIsLessThan<double>(double x, double y);
  template bool NormIsLessThan<CD>(CD x, double y);
  template int ErfcAtN<double>(double x, unsigned int n,
			       double& y, double& yn);
  template int ErfcAtN<CD>(CD x, unsigned int n, CD& y, CD& yn);
  template int Erfc<double>(double, double&, ErfcCalcData&);
  template int Erfc<CD>(CD, CD&, ErfcCalcData&);
  template int Exp2ErfcAtN<double>(double x, unsigned int n,
				   double&, double&);
  template int Exp2ErfcAtN<CD>(CD x, unsigned int n,
			       CD&, CD&);  

  template int Exp2Erfc<double>(double x, double& y,
				ErfcCalcData& data);
  template int Exp2Erfc<CD>(CD x, CD& y, ErfcCalcData& data);  
			    
}
