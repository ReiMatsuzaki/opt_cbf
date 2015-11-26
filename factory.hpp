#ifndef FACTORY_HPP
#define FACTORY_HPP

#include <complex>
#include <vector>
#include <stdexcept>
#include <Eigen/Core>

/**
 * this factory object create IOptTarget, IOptimizer, LinearComb<Prim> object
 * from KeyValues object.
 */

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
  template<class F> class IOptimizer;
  class IOptTarget;

  class IFactory {
  protected:
    KeysValues* kv_;
  public:
    explicit IFactory(const KeysValues&);
    virtual ~IFactory();
    IOptTarget* OptTarget() const; 
    virtual IOptTarget* optTarget() const = 0;
    virtual IOptimizer<CD>* Optimizer() const = 0;
    /**
     * getter for orbital exponents of basis function
     */
    virtual void GetZs(VectorXcd* zs) const = 0;
    /**
     * getter for the number of basis
     */
    int BasisSize() const;

  };

  // ============== Mono ================================
  template<class PrimBasis, class PrimDriv>
  class FactoryMono : public IFactory {
  public:
    FactoryMono(const KeysValues&);
    ~FactoryMono();
    vector<PrimBasis>* BasisSet() const;
    IOptTarget* optTarget() const;
    IOptimizer<CD>* Optimizer() const; 
    void GetZs(VectorXcd* zs) const;
  };

  // ============== EvenTempered ========================
  template<class PrimBasis, class PrimDriv>
  class FactoryEvenTemp : public IFactory {
  public:
    FactoryEvenTemp(const KeysValues&);
    ~FactoryEvenTemp();
    vector<PrimBasis>* BasisSet() const;
    IOptTarget* optTarget() const;
    IOptimizer<CD>* Optimizer() const;
    void GetZs(VectorXcd* zs) const;
  };

}

#endif
