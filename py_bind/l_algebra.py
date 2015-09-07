import numpy as np
from numpy import dot, outer
import scipy.linalg as la
"""
Same functions with ../l_algebra, but implemented in scipy.linalg
"""

def calc_a_Aj_b(a, A10, b):
    A10_b = dot(A10, b)
    A10_a = dot(A10, a)
    return a * A10_b + b * A10_a

def calc_a_Ai_B_Aj_b(a, A10, B, b):
    t1 = outer(a, A10.dot(b)) * dot(A10, B)
    t2 = outer(dot(A10, a), dot(A10, b)) * B
    t3 = dot(dot(A10, B), A10.T) * outer(a, b)
    t4 = outer(dot(A10, a), b) * dot(B, A10.T)
    return t1 + t2 + t3 + t4

def calc_a_Aij_a(a, A20, A11):
    t1 = outer(a, a) * A11
    t2 = np.diag(a * dot(A20, a))
    return 2.0 * (t1 + t2)

def calc_ai_A_Bj_b(a1, A, B10, b):
    t1 = outer(a1, dot(B10, b)) * A
    t2 = outer(a1, b) * dot(A, B10.T)
    return t1 + t2

def calc_ai_b(a1, b):
    return a1 * b

def calc_val_grad_hess_aAm1a(a0, a1, a2, A00, A10, A20, A11):
    """
    a_i = (u_i(z_i), a)
    A_ij= (u_i(z_i), A u_j(z_j))
    compute value, gradient and hessian of aA^{-1}a
    """

    # A00^{-1}a
    Am1a = la.solve(A00, a0)
    
    # value
    val = dot(a0, Am1a)

    if a1 == None or a2 == None or A10 == None or A20 == None or A11 == None:
        return (val, None, None)

    # gradient
    t1 = -1 * calc_a_Aj_b(Am1a, A10, Am1a)
    t2 = +2 * calc_ai_b(a1, Am1a)
    grad = t1 + t2

    # hessian
    A_inv = la.inv(A00)
    t1  = np.diag(2 * a2 * Am1a)
    t2  = +2 * outer(a1, a1) * A_inv
    t3  = -1 * calc_a_Aij_a(Am1a, A20, A11)
    t4  = -2 * calc_ai_A_Bj_b(a1, A_inv, A10, Am1a)
    t5  = t4.T
    t6  = 2 * calc_a_Ai_B_Aj_b(Am1a, A10, A_inv, Am1a)
    hess = t1 + t2 + t3 + t4 + t5 + t6

    return (val, grad, hess)


