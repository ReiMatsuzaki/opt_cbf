import unittest
import l_algebra
import numpy as np
from numpy import dot, outer
import sys
sys.path.append('../../l2func')
import l2func as l2
import opt_cbf

class TestCalculations(unittest.TestCase):
    def setUp(self):
        pass

    def test_a_Aj_b(self):
        # copied from ../utest.cpp
        """
        // A^0 = (2.0, 2.0)
        //       (2.0, 0.0)
        // A^1 = (0.0,   1.5)
        //       (1.5,   0.2)
        // aA^0b = (0.1,0.2) (1.6, 0.6) = 0.16 + 0.12 = 0.28
        // aA^1b = (0.1,0.2) (0.75, 0.45+0.1) = 0.075+0.11=0.185
        MatrixXd A10(2,2); A10 << 1.0, 2.0, 1.5, 0.1;
        VectorXd a(2), b(2);  a << 0.1, 0.2;  b << 0.3, 0.5;
        VectorXd res(2);
        Calc_a_Aj_b(a, A10, b, &res);

        EXPECT_DOUBLE_EQ(0.28, res(0));
        EXPECT_DOUBLE_EQ(0.185, res(1));        
        """
        mat_A10 = np.array([[1.0, 2.0], [1.5, 0.1]])
        vec_a  = np.array([0.1, 0.2])
        vec_b  = np.array([0.3, 0.5])
        res = l_algebra.calc_a_Aj_b(vec_a, mat_A10, vec_b)
        self.assertAlmostEqual(0.28,  res[0])
        self.assertAlmostEqual(0.185, res[1])
        
    def test_a_Ai_B_Aj_b(self):
        # copied from ../utest.cpp
        """
        MatrixXd A10(2,2);   A10 << 1.002, 1.008, 1.003, 1.016;
        VectorXd a(2), b(2); a   << 0.2, 0.4; b << 0.3, 0.1;
        MatrixXd B(2,2);     B   << 1.0, 0.1, 0.1, 0.3;

        MatrixXd res(2,2);
        Calc_a_Ai_B_Aj_b(a, A10, B, b, &res);

        MatrixXd A_d0_d0(2, 2); A_d0_d0 << 2.004, 1.008, 1.008, 0.0;
        EXPECT_DOUBLE_EQ( a.transpose() * A_d0_d0 * B * A_d0_d0 * b,res(0, 0));
        """
        A10 = np.array( [[1.002, 1.008], [1.003, 1.016]])
        a   = np.array( [0.2, 0.4])
        b   = np.array( [0.3, 0.1])
        B   = np.array( [[1.0, 0.1], [0.1, 0.3]])
        res = l_algebra.calc_a_Ai_B_Aj_b(a, A10, B, b)
        A_d0_d0 = np.array([[2.004, 1.008], [1.008, 0.0]])
        self.assertAlmostEqual( dot(a, dot(dot(dot(A_d0_d0, B), A_d0_d0), b)), res[0][0])

    def test_a_Aij_b(self):
        # copied from ../utest.cpp
        """
        MatrixXd A20(2,2); A20 << 6.6, 6.6, 7.2, 7.2;
        MatrixXd A11(2,2); A11 << 1.0, 1.0, 1.0, 1.0;
        VectorXd a(2)    ; a   << 0.2, 0.4;

        MatrixXd res(2,2);
        Calc_a_Aij_a(a, A20, A11, &res);
        EXPECT_DOUBLE_EQ(1.664, res(0,0));
        """
        A20 = np.array([[ 6.6, 6.6], [7.2, 7.2]])
        A11 = np.array([[1.0, 1.0],  [1.0, 1.0]])
        a   = np.array([ 0.2, 0.4])
        res = l_algebra.calc_a_Aij_a(a, A20, A11)
        self.assertAlmostEqual(1.664, res[0][0])

    def test_grad(self):
        us   = [ l2.STO(1.0, 2, z) for z in [1.1, 1.2]]
        driv = l2.STOs()
        driv.add_one(1.0, l2.STO(1.0 / np.sqrt(2.0), 2, 1.0))

        h_atom = l2.h_like_atom('2p')
        l_op = h_atom.h_minus_energy_sto(0.5)

        (val, grad, hess, datas) = opt_cbf.val_grad_hess(us, driv, l_op)

        """
        These values are taken from 2015/9/1_**.ipynb
        """

        self.assertAlmostEqual(0.608906838305313,  datas['a0'][0])
        self.assertAlmostEqual(0.599798428001011,  datas['a0'][1])

        self.assertAlmostEqual(-0.445,             datas['A00'][0][0])
        self.assertAlmostEqual(-0.413041532083526, datas['A00'][1][0])
        self.assertAlmostEqual(-0.413041532083526, datas['A00'][0][1])
        self.assertAlmostEqual(-0.380000000000000, datas['A00'][1][1])

        self.assertAlmostEqual(-0.478372504129766, val)
        self.assertAlmostEqual(-1.01381241031413, grad[0])
        self.assertAlmostEqual(-0.915041952889624, grad[1])

        self.assertAlmostEqual(6.42894956208815, hess[0][0])
        self.assertAlmostEqual(5.10430865981922, hess[0][1])
        self.assertAlmostEqual(5.10430865981922, hess[1][0])
        self.assertAlmostEqual(5.77073843522357, hess[1][1])


class TestOpt(unittest.TestCase):
    def test_1basis_1skp(self):
        us = [ l2.STO(1.0, 2, 0.5-0.5j) ]
        res = opt_cbf.optimize_hydrogen_pi(us, '1s->kp', 0.9)
        self.assertTrue(res[0])
        zs = [ basis.z for basis in res[-1]]
        self.assertAlmostEqual(1.1117640506-0.3673558953j, zs[0])
        
    def test_3basis_1skp(self):
        us = [ l2.STO(1.0, 2, z) for z 
               in [0.96-0.008j, 1.04-0.71j, 0.45-1.37j]]
        res = opt_cbf.optimize_hydrogen_pi(us, '1s->kp', 0.9, eps = 0.000001, max_iter = 100)
        self.assertTrue(res[0])
        zs = [ basis.z for basis in res[-1]]
        self.assertAlmostEqual(0.9899070351-0.0123382063j, zs[0])
        self.assertAlmostEqual(1.0418217968-0.7127594141j, zs[1])
        self.assertAlmostEqual(0.4509924544-1.3744734865j, zs[2])

        (psi, cs, a) = opt_cbf.solve(res[-1], channel = '1s->kp', energy = 0.9)
        self.assertAlmostEqual(0.7825445397-0.1236019852j, cs[0] / (-np.sqrt(3.0)))
        self.assertAlmostEqual(0.2031470249-0.0913057736j, cs[1] / (-np.sqrt(3.0)))
#        self.assertAlmostEqual(0.0460737405+0.0097070218j, cs[2] / (-np.sqrt(3.0)))
        self.assertAlmostEqual(0.7825445397-0.1236019852j, psi.coef_i(0) / (-np.sqrt(3.0)))
        self.assertAlmostEqual(0.2031470249-0.0913057736j, psi.coef_i(1) / (-np.sqrt(3.0)))
        #        self.assertAlmostEqual(0.0460737405+0.0097070218j, cs[2] / (-np.sqrt(3.0)))
    def test_2gto_1skp(self):
        # us = [ l2.STO(2, z) for z in [(0.0361962-0.0271314j), (0.148847-0.145461j)] ]
        us = [ l2.GTO(1.0, 2, z) for z in [(0.036-0.027j), (0.15-0.14j)] ]
        eps = 0.0000001
        opt_res = opt_cbf.optimize_hydrogen_pi(us, '1s->kp', 0.5, eps = eps)
        self.assertTrue(opt_res[-1])
        zs = [ basis.z for basis in opt_res[-1]]
        self.assertTrue( abs(0.0361962-0.0271314j - zs[0]) < eps * 10)
        self.assertTrue( abs(0.148847-0.145461j   - zs[1]) < eps * 10)

        """
    def test_solve_by_cbf(self):
        us = [ l2.STO(1.0, 2, z) for z in [1.0, 0.6-0.4j]]
        (psi, cs, a) = opt_cbf.solve(us, channel = '1s->kp', energy = 0.5)
        print psi
"""

if __name__ == '__main__':
    unittest.main()
