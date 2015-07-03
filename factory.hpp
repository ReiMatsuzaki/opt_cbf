#ifndef FACTORY_HPP
#define FACTORY_HPP

#include <complex>
#include <vector>
#include <stdexcept>
#include <Eigen/Core>

namespace {
  using std::vector;
  using std::string;
  using Eigen::VectorXcd;
}

class KeysValues;

namespace l2func {
  template<class F, int m> class ExpBasis;
  typedef ExpBasis<std::complex<double>, 1> CSTO;
  typedef ExpBasis<std::complex<double>, 2> CGTO;
  
}

namespace opt_cbf_h {

  // =============== Create Factory ====================
  class IFactory;
  IFactory* CreateFactory(const KeysValues&);

  // ============== Exceptions =========================
  class InvalidBasis : public std::runtime_error {
  private:
    string basis_type_;
  public:
    InvalidBasis(string msg, string basis_type);
    ~InvalidBasis() throw();
    const char* what() const throw();
  };
  class InvalidDriv : public std::runtime_error {
  private:
    string ch_, di_;
  public:
    InvalidDriv(string msg, string ch, string di);
    ~InvalidDriv() throw();
    const char* what() const throw();
  };

  // =============== Interface =========================
  template<class Prim> class HAtomPI;
  template<class F> class IOptimizer;
  class IOptTarget;

  class IFactory {
  protected:
    KeysValues* kv_;
  public:
    explicit IFactory(const KeysValues&);
    virtual ~IFactory();
    virtual vector<l2func::CSTO>* STOSet() const = 0;
    virtual vector<l2func::CGTO>* GTOSet() const = 0;
    HAtomPI<l2func::CSTO>* HAtomPiSTO() const;
    HAtomPI<l2func::CGTO>* HAtomPiGTO() const;
    IOptTarget* OptTarget() const; 
    void SetZs(VectorXcd* zs) const;
    virtual IOptimizer<CD>* Optimizer() const = 0;
  };

  // ============== Mono ================================
  class FactoryMono : public IFactory {
  public:
    FactoryMono(const KeysValues&);
    ~FactoryMono();
    vector<l2func::CSTO>* STOSet() const;
    vector<l2func::CGTO>* GTOSet() const;
    IOptimizer<CD>* Optimizer() const;    
  };

  // ============== EvenTempered ========================
  class FactoryEvenTemp : public IFactory {
  public:
    FactoryEvenTemp(const KeysValues&);
    ~FactoryEvenTemp();
    vector<l2func::CSTO>* STOSet() const;
    vector<l2func::CGTO>* GTOSet() const;
    IOptimizer<CD>* Optimizer() const;
  };

}

#endif
