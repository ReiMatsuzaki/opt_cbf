#include "erfc.hpp"
#include "fact.hpp"
#include "prim.hpp"

namespace {
  using erfc_mori::ErfcCalcData;

}

namespace l2func {

  // =========== basic calc ===================
  const CD operator*(const CD& a, int b) {
    return CD(a.real() * b, a.imag() * b); }
  const CD operator*(int b, const CD& a) {
    return CD(a.real() * b, a.imag() * b); }  
  int int_pow(int base, unsigned int expo) {

    int acc = 1;
    for(unsigned int i = 0; i < expo; i++) 
      acc *= base;
    return(acc);
    
  }
  template<class F> F STO_Int(F z, int n) {
    return pow(z, -n-1.0) * (1.0 * fact::Factorial(n));
  }
  template<class F> F GTO_Int(F z, int n) {
    
    F res;
    if(n % 2 == 0) {
      
      int nn = n/2;
      res = fact::DoubleFactorial(2*nn-1) * sqrt(M_PI) /
	(F(int_pow(2, nn+1)) * pow(sqrt(z), 2*nn+1));
      
    } else {
      
      int nn = (n-1)/2;
      res = F(fact::Factorial(nn)) / (F(2) * pow(z, nn+1));
      
    }
    return res;
  }  

  // ============ ExpBasis ===================
  template<class F, int m>
  ExpBasis<F,m>::ExpBasis() : c_(0), n_(0), z_(F(1.0)) {}
  template<class F, int m>
  ExpBasis<F,m>::ExpBasis(int _n, F _z):c_(1), n_(_n),z_(_z) {}
  template<class F, int m>
  ExpBasis<F,m>::ExpBasis(F _c, int _n, F _z) :
    c_(_c), n_(_n), z_(_z) {}
  template<class F, int m>
  ExpBasis<F,m>::ExpBasis(int _n, F _z, ENormalized):
    c_(1), n_(_n), z_(_z) {

      F c2;
      if(m == 1)
	c2 = STO_Int<F>(2 * z_, 2 * n_);
      else
	c2 = GTO_Int<F>(2 * z_, 2 * n_);
      
      c_ = F(1) / sqrt(c2);
    }
  
  template<class F, int m>
  std::ostream& operator << (std::ostream& os, 
			     const ExpBasis<F,m>& a) {
    //    int m = Prim::exp_power;
    //    typename Prim::Field c = a.c();
    //typename Prim::Field z = a.z();
    F c = a.c();

    if(m == 1) 
      os << (char *)"STO(" ;
    else
      os << (char *)"GTO(";
    char* comma = (char*)", ";
    os << c;
    os << comma << a.n() << comma << a.z() << (char*)")";
    return os;
  }
  template std::ostream& operator <<<CD,1> (std::ostream& os, const CSTO& a);
  template std::ostream& operator <<<double,1> (std::ostream& os, const RSTO& a);
  template std::ostream& operator <<<CD,2> (std::ostream& os, const CGTO& a);
  template std::ostream& operator <<<double,2> (std::ostream& os, const RGTO& a);


  // =========== inner product ===================
  template<class F> F sto_gto_int_0(F as, F ag) {
    F erfcVal, expVal, sqrtPi,pi,res;
    ErfcCalcData data;
    Erfc(as/(2*sqrt(ag)),erfcVal,data);
    expVal=exp(as*as/(4*ag));
    pi=M_PI;
    sqrtPi=sqrt(pi);
    res = (erfcVal*expVal*sqrtPi)/(2*sqrt(ag));
    return (res);
  }
  template<class F> F sto_gto_int_1(F as, F ag) {
    F exp2erfc, sqrtPi, sqrt_ag, pi, res;
    ErfcCalcData data;
    sqrt_ag = sqrt(ag);
    pi = M_PI;
    sqrtPi = sqrt(pi);
    Exp2Erfc(as/(2*sqrt_ag), exp2erfc, data);

    res = F(1)/(2*ag) - (as * exp2erfc * sqrtPi) /
      (4*ag*sqrt_ag);  
  
    return(res);
  }
  template<class F> F sto_gto_int_2(F as, F ag) {
    F erfcVal, expVal, sqrtPi,pi,res;
    ErfcCalcData data;
    Erfc(as/(2*sqrt(ag)),erfcVal,data);
    expVal=exp(as*as/(4*ag));
    pi=M_PI;
    sqrtPi=sqrt(pi);
    res = (-2*sqrt(ag)*as + (2*ag + pow(as,2))*erfcVal*expVal*sqrtPi)/(8*pow(sqrt(ag),5));

    return (res);
  }
  template<class F> F sto_gto_int_3(F as, F ag) {
    F sqrtPi,pi,res, sqrt_ag, exp2erfcVal;
    ErfcCalcData data;

    sqrt_ag = sqrt(ag);
    pi=M_PI; 
    sqrtPi=sqrt(pi);
    Exp2Erfc(as/(2*sqrt_ag), exp2erfcVal, data);
    
    res = (4*ag + pow(as,2))/(8*pow(ag,3)) 
      -(as*(6*ag + pow(as,2)) * exp2erfcVal * sqrtPi) / 
      (16 * sqrt_ag * ag * ag * ag);

    return (res);
  }
  template<class F> F sto_gto_int_4(F as, F ag) {
    F sqrtPi,pi,res, sqrt_ag, exp2erfcVal;
    ErfcCalcData data;

    sqrt_ag = sqrt(ag);
    pi=M_PI; 
    sqrtPi=sqrt(pi);
    Exp2Erfc(as/(2*sqrt_ag), exp2erfcVal, data);
    
    res = (-2*sqrt_ag*as*(10*ag + as*as)
	   + (12*ag*ag + 12*ag*as*as + as*as*as*as)*exp2erfcVal*sqrtPi)/(32*ag*ag*ag*ag*sqrt_ag);
    
    return (res);
  }
  template<class F> F sto_gto_int_5(F as, F ag) {
    F sqrtPi,pi,res, sqrt_ag, exp2erfcVal;
    ErfcCalcData data;

    Exp2Erfc(as/(2*sqrt(ag)), exp2erfcVal, data);
    pi=M_PI;
    sqrtPi=sqrt(pi);
    sqrt_ag = sqrt(ag);
    F x = sqrt_ag;
 
    res = (2*sqrt_ag*(2*ag + as*as)*(16*ag + as*as) -
	   as*(60*ag*ag + 20*ag*as*as + as*as*as*as)*
	   exp2erfcVal*sqrtPi)/
      (64*x*ag*ag*ag*ag*ag);
    return (res);
  }
  template<class F> F sto_gto_int_6(F as, F ag) {

    F sqrtPi,pi,res, sqrt_ag, exp2erfcVal;
    ErfcCalcData data;

    sqrt_ag = sqrt(ag);
    pi=M_PI; 
    sqrtPi=sqrt(pi);
    F a = as/(2*sqrt_ag);
    Exp2Erfc(a, exp2erfcVal, data);

    if(! data.convergence) {
      string msg("is not convergence in sto_gto_int_6");
      throw msg;
    }
    
    res = (-2*sqrt_ag*as*(6*ag + as*as)*(22*ag + as*as) +
	   (120*ag*ag*ag + 180*ag*ag*as*as + 30*ag*pow(as,4) + 
	  pow(as,6))*exp2erfcVal*sqrtPi)/
      (128*sqrt_ag*ag*ag*ag*ag*ag*ag);

    return (res);
  }
  template<class F> F sto_gto_int_7(F as, F ag) {
    
    F sqrtPi,pi,res, sqrt_ag, exp2erfcVal;
    ErfcCalcData data;

    sqrt_ag = sqrt(ag);
    pi=M_PI; 
    sqrtPi=sqrt(pi);
    Exp2Erfc(as/(2*sqrt_ag), exp2erfcVal, data);
  
    res = (2*sqrt(ag)*(384*pow(ag,3) + 348*pow(ag,2)*pow(as,2) + 40*ag*pow(as,4) + pow(as,6)) - as*(840*pow(ag,3) + 420*pow(ag,2)*pow(as,2) + 42*ag*pow(as,4) + pow(as,6))*exp2erfcVal*sqrtPi)/
    (256*sqrt_ag * ag* ag* ag* ag* ag* ag* ag);

    return (res);
  }
  template<class F> F sto_gto_int_8(F as, F ag) {
    F sqrtPi,pi,res, sqrt_ag, exp2erfcVal;
    ErfcCalcData data;

    sqrt_ag = sqrt(ag);
    pi=M_PI; 
    sqrtPi=sqrt(pi);
    Exp2Erfc(as/(2*sqrt_ag), exp2erfcVal, data);
    
    res = (-2*sqrt(ag)*as*(2232*pow(ag,3) + 740*pow(ag,2)*pow(as,2) + 54*ag*pow(as,4) + pow(as,6)) + (1680*pow(ag,4) + 3360*pow(ag,3)*pow(as,2) + 840*pow(ag,2)*pow(as,4) + 56*ag*pow(as,6) + pow(as,8))*exp2erfcVal*sqrtPi)/(512*sqrt_ag*pow(ag,8));

    return (res);
  }
  template<class F> F sto_gto_int_9(F as, F ag) {

    F sqrtPi,pi,res, sqrt_ag, exp2erfcVal;
    ErfcCalcData data;

    sqrt_ag = sqrt(ag);
    pi=M_PI; 
    sqrtPi=sqrt(pi);
    Exp2Erfc(as/(2*sqrt_ag), exp2erfcVal, data);

    res = (2*sqrt_ag*(6144*pow(ag,4) + 7800*pow(ag,3)*pow(as,2) + 1380*pow(ag,2)*pow(as,4) + 70*ag*pow(as,6) + pow(as,8)) - as*(15120*pow(ag,4) + 10080*pow(ag,3)*pow(as,2) + 1512*pow(ag,2)*pow(as,4) + 72*ag*pow(as,6) + pow(as,8))*exp2erfcVal*sqrtPi)/(1024*sqrt_ag*pow(ag,9));


    return (res);
  }
  template<class F> F STO_GTO_Int(F as, F ag, int n) {
    F res;
    switch(n) {
    case 0:
      res = sto_gto_int_0(as, ag);
      break;
    case 1:
      res = sto_gto_int_1(as, ag);
      break;
    case 2:
      res = sto_gto_int_2(as, ag);
      break;
    case 3:
      res = sto_gto_int_3(as, ag);
      break;
    case 4:
      res = sto_gto_int_4(as, ag);
      break;
    case 5:
      res = sto_gto_int_5(as, ag);
      break;
    case 6:
      res = sto_gto_int_6(as, ag);
      break;
    case 7:
      res = sto_gto_int_7(as, ag);
      break;
    case 8:
      res = sto_gto_int_8(as, ag);
      break;
    case 9:
      res = sto_gto_int_9(as, ag);
      break;      
    default:
      string msg;
      msg = "this is not supported in STO_GTO_int";
      throw msg;
    }

    return res;
  }
  
  template<class F>
  F CIP(const ExpBasis<F, 1>& a, const ExpBasis<F, 1>& b) {
    return STO_Int(a.z() + b.z(), a.n() + b.n()) * a.c() * b.c();
  }
  template<class F>
  F CIP(const ExpBasis<F, 1>& a, const ExpBasis<F, 2>& b) {
    return STO_GTO_Int(a.z(), b.z(), a.n() + b.n()) * a.c() * b.c();
  }
  template<class F>  
  F CIP(const ExpBasis<F, 2>& a, const ExpBasis<F, 1>& b) {
    return CIP(b, a);
  }
  template<class F>
  F CIP(const ExpBasis<F, 2>& a, const ExpBasis<F, 2>& b) {
    return GTO_Int(a.z() + b.z(), a.n() +
		   b.n()) * a.c() * b.c();
  }  

  
  // explicit instance
  typedef std::complex<double> CD;
  
  template double STO_Int<double>(double z, int n);
  template CD STO_Int<CD>(CD z, int n);
  template double GTO_Int<double>(double z, int n);
  template CD GTO_Int<CD>(CD z, int n);  
  
  template class ExpBasis<double, 1>;
  template class ExpBasis<double, 2>;
  template class ExpBasis<CD, 1>;
  template class ExpBasis<CD, 2>;

  template double CIP(const RSTO&, const RSTO&);
  template double CIP(const RGTO&, const RSTO&);
  template double CIP(const RSTO&, const RGTO&);
  template double CIP(const RGTO&, const RGTO&);
  template CD CIP(const CSTO&, const CSTO&);
  template CD CIP(const CGTO&, const CSTO&);
  template CD CIP(const CSTO&, const CGTO&);
  template CD CIP(const CGTO&, const CGTO&);
}
