import unittest
import nnewton
import numpy as np

def func(k, xs):
    x = xs[0]
    return np.cos(x-k)

class TestNNewton(unittest.TestCase):
    def setUp(self):
        pass

    def test_nnewton(self):
        k0 = 0.2
        (convq, xs_o, val, grad) = nnewton.nnewton(lambda xs:func(k0, xs), [0.0])
        self.assertTrue(convq)
        self.assertAlmostEqual(k0, xs_o[0])

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
