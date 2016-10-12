#ifndef RUN_TEMPLATE_TEMPLATE_H
#define RUN_TEMPLATE_TEMPLATE_H

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

namespace opt_cbf_h {

  template<class R, class B, class S, class Op>
  struct Run {
    void Parse(int argc, char *argv[]) {

      // parse
      cmdline::parser parser;
      parser.add<vector<cd> >("zs", 'z', "complex", true);
      parser.add<double>("w", 'w', "photon energy", true);
      parser.add<int>("maxit", 'm', "max iteration", false, 100);
      parser.add<double>("eps", 'e', "convergence epsilon", false, 0.0000001);
      parser.parse_check(argc, argv);

      // extract values
      vector<cd> zs = parser.get<vector<cd> >("zs");
      double w      = parser.get<double>("w");
      double ene0   = -0.5;
      int   maxit   = parser.get<int>("maxit");
      double eps    = parser.get<double>("eps");

      // ---- print value ----
      cout << "Input" << endl;
      int i(0);
      for(vector<cd>::const_iterator it = zs.begin(); it != zs.end(); ++it) {
	cout << "in_zs_" << i << ": " << *it << endl;
	i++;
      }
      cout << "in_w:"   << w  << endl;
      cout << "in_ene0: " << ene0 << endl;
      cout << "in_maxit: "   << maxit << endl;
      cout << "in_eps: "     << eps   << endl;
      cout << endl;

      // create object
      zs0 = Eigen::Map<VectorXcd>(&zs[0], zs.size());

      HLikeAtom<cd> hatom;
      HLength<1, 0, 1, cd> mu_phi(hatom);
      HminusEOp<1, cd> lop(hatom, ene0 + w);

    }

    
    
    VectorXcd zs0;
    OptTarget<R,B,S,Op>* opt_target;
    IOptimizer<cd>* optimizer;
    OptRes<cd> res;

  };
}


#endif
