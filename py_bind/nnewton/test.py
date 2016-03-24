import unittest
from nnewton import *
import numpy as np

def func(k, xs):
    x = xs[0]
    return np.cos(x-k)

def func_grad_hess(k, xs):
    x = xs[0]
    a = 1.5
    return (np.cos(a*x-k),
            np.array([-a*np.sin(a*x-k)]),
            np.array([[-a*a*np.cos(a*x-k)]]))

class TestNNewton(unittest.TestCase):
    def setUp(self):
        pass

    def test_pd(self):
        x = [1.1]
        (val, grad, hess) = func_grad_hess(0.3, x)
        h = 0.001
        self.assertAlmostEqual(grad, num_pd(func, x, h, 1, 0, "r1")+0.1)
        self.assertAlmostEqual(grad, num_pd(func, x, h, 1, 0, "c1"))

    def test_pd_(self):
        x = [1.1]
        (val, grad, hess) = func_grad_hess(0.3, x)
        h = 0.001
        self.assertAlmostEqual(grad, num_pd(func, x, h, 1, 0, "r1")+0.1)
        self.assertAlmostEqual(grad, num_pd(func, x, h, 1, 0, "c1"))        

    def test_nnewton(self):
        k0 = 0.2
        (convq, xs_o, val, grad) = nnewton(lambda xs:func(k0, xs), [0.0])
        self.assertTrue(convq)
        self.assertAlmostEqual(k0, xs_o[0])

    def test_nnewton_analytic(self):
        k0 = 0.2
        target = lambda xs:func_grad_hess(k0, xs)

        self.assertTrue(nnewton.is_reasonable_val_grad_hess(target, [0.1], method="c1", h=0.001))

        (convq, xs_o, val, grad) = nnewton(target, [0.1], func="vgh", show_lvl=1)
        self.assertTrue(convq)
        self.assertAlmostEqual(k0/1.5, xs_o[0])

    def test_nnewton_params(self):
        ks = [-0.1, 0.0, 0.1]
        res = nnewton.nnewton_params(func, ks, (0.0, [0.03]))
        (k, xs_o, val) = res[0]
        self.assertAlmostEqual(k, -0.1)
        self.assertAlmostEqual(k, xs_o[0])
        (k, xs_o, val) = res[1]
        self.assertAlmostEqual(k, 0.0)
        self.assertAlmostEqual(k, xs_o[0])
        (k, xs_o, val) = res[2]
        self.assertAlmostEqual(k, 0.1)
        self.assertAlmostEqual(k, xs_o[0])


if __name__ == '__main__':
    unittest.main()
