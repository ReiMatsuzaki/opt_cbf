#ifndef FACT_HPP
#define FACT_HPP

namespace fact {
  
  int Factorial(int n) {

    if(n < 0) {
      throw std::invalid_argument("n must be zero or positive in Factorial.");
    }

    int acc = 1;
    for(int i = n; i > 0; i--)
      acc *= i;

    return acc;
  }

  int DoubleFactorial(int n) {

    if(n < 0) {
      throw std::invalid_argument("n must be zero or positive in DoubleFactorial.");
    }


    int acc = 1;
    for(int i = n; i > 0; i-=2) {
      acc *= i;
    }
    
    return acc;
  }

  /*

  template<int n>
  struct FactStruct {
    enum { value = n * FactStruct<N-1>::value };
  };

  template<>
  struct FactStruct<0> {
    enum { value = 1 };
  };

  class FactData {
  private:
    vector<int> fact_list_;
    FactData(int num) {
      FactData.resize(num);
      for(int i = 0 ; i < num; i++)
    }
  };

  int Fact(int n) {
    
  }
  */
}

#endif
