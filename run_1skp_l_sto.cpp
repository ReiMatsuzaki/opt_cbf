#include <vector>

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <Eigen/Core>

#include "opt_target.hpp"
#include "opt_target_impl.hpp"
#include "opt.hpp"

#include "external/cmdline.h"

using namespace std;
using namespace opt_cbf_h;
using namespace l2func;
typedef complex<double> cd;

int main(int argc, char *argv[]) {

  // ---- options parser ----
  cmdline::parser parser;
  parser.add<vector<cd> >("zs", 'z', "complex", true);
  parser.add<double>("w", 'w', "photon energy", true);
  parser.add<int>("maxit", 'm', "max iteration", false, 100);
  parser.add<double>("eps", 'e', "convergence epsilon", false, 0.0000001);
  parser.parse_check(argc, argv);

  // ---- set values ----
  vector<cd> zs = parser.get<vector<cd> >("zs");
  double w      = parser.get<double>("w");
  double ene0   = -0.5;
  int   maxit   = parser.get<int>("maxit");
  double eps    = parser.get<double>("eps");

  // ---- print value ----
  int i(0);
  for(vector<cd>::const_iterator it = zs.begin(); it != zs.end(); ++it) {
    cout << "in_zs_" << i << ": " << *it << endl;
    i++;
  }
  cout << "in_w:"   << w  << endl;
  cout << "in_ene0: " << ene0 << endl;
  cout << "in_maxit: "   << maxit << endl;
  cout << "in_eps: "     << eps   << endl;

  // ---- build calculation objects ----
  VectorXcd zs0 = Eigen::Map<VectorXcd>(&zs[0], zs.size());

  HLikeAtom<cd> hatom;
  HLength<1, 0, 1, cd> mu_phi(hatom);
  HminusEOp<1, cd> lop(hatom, ene0 + w);

  vector<NormalCSTO> basis_set;
  for(vector<cd>::const_iterator it = zs.begin(); it != zs.end(); ++it) 
    basis_set.push_back(NormalCSTO(2, *it));
  

  // ---- optimization ----
  OptSTOLength *opt_target = 
    new OptSTOLength(mu_phi.value, basis_set, mu_phi.value, lop.value);
  IOptimizer<cd> *opt = new OptimizerNewton<cd>(maxit, eps);
  OptRes<cd> opt_res = opt->Optimize(bind(&IOptTarget::Compute, opt_target,
					  _1, _2, _3, _4), zs0);


  // ---- output ----
  cout << "convergence: " << opt_res.convergence << endl;
  cout << opt_res.z << endl;

  delete opt_target;
  delete opt;

  return 0;
}
