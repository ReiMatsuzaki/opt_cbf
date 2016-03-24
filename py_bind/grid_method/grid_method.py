import os
import sys
import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg

sys.path.append("../../l2func/py_bind")
from l2func import *

# ==== utility ====
def flatten(xss):
    return reduce(lambda a,b:a+b, xss)

def laplacian_mat(n):
    data = [1, -2, 1]*n
    i = flatten([[k,k,k] for k in range(n)])
    j = flatten([[k-1, k, k+1] for k in range(n)])
    return scipy.sparse.coo_matrix((data[1:-1], (i[1:-1], j[1:-1])))

def bc_outgoing_mat(n, h, k):
    d = [1.0, 2.0j*k*h]
    i = [n-1, n-1]
    j = [n-2, n-1]
    return scipy.sparse.coo_matrix((d, (i, j)))

# ==== solver ====
def solve_driv(v, ene, s, n, h):
    """
    solve (T+V-E)f=s with outgoing boundary condition
    
    Inputs
    ------
    v: lambda scalar->scalar
      potential
    ene: scalar
      E in (T+V-E)
    s: lambda scalar->scalar
      driven term
    n: integer 
      number of grid points
    h: real
      grid width
    """

    xs = np.array([(k+1)*h for k in range(n)])
    h2 = h*h
    k = np.sqrt(2.0*ene)
    
    vs = [v(x)-ene for x in xs]

    mat = laplacian_mat(n) -2.0 * h2 * scipy.sparse.diags(vs, 0) + bc_outgoing_mat(n, h, k)
    vec = np.array([-2.0*h*h*s(x) for x in xs])

    ys = scipy.sparse.linalg.spsolve(mat, vec)
    return (xs, ys)

def solve_h_length(ene, channel, n, h):
    """
    solve (H-E)f=s with outgoing boundary condition.

    Inputs
    ------
    ene : scalar
       E in operator (H-E)
    channel : integer tuple 
       (n0, l0, l1) where n0 is initial quantum number, l0 and l1 are 
       initial and final angular quantum number 
    n : number of grids
    h : grid width

    Returns
    -------
    xs : array
    ys : solution
    """

    (n0, l0, l1) = channel
    v = lambda x: -1.0/x + l1*(l1+1)*0.5/(x*x)
    driv_term = HAtom(1.0).length(n0, l0, l1)
    print "driv_term>>>"
    print driv_term
    print "driv_term<<<"
    s = lambda x: driv_term.at(x)
    return solve_driv(v, ene, s, n, h)


