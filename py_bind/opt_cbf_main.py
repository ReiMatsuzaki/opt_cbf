import sys
import copy
import cmath
import numpy as np
from numpy import sqrt
from scipy.special import gamma
import scipy.linalg as la
import l_algebra
import l2func as l2
from utils import *

# ======= Utils =========
def get_sym_ip(type1, type2):
    if   type1 == l2.GTO and type2 == l2.GTO:
        return l2.sym_ip_gg
    elif type1 == l2.GTO and type2 == l2.STO:
        return l2.sym_ip_gs
    elif type1 == l2.STO and type2 == l2.GTO:
        return l2.sym_ip_sg
    elif type1 == l2.STO and type2 == l2.STO:
        return l2.sym_ip_ss
    else:
        print 'unsupported type! in opt_cbf.py::val_grad_hess'
        sys.exit()

def convergence_q(eps):
    def is_small(xs):
        return max([abs(x) for x in xs]) < eps
    def __func__(grad, dz):
        return is_small(grad) and is_small(dz)
    return __func__
        
def channel_to_op_and_driv(channel, ene):

    [init, final] = channel.split("->")    
    if(channel   == "1s->kp"):
        (n0, l0, l1) = (1, 0, 1)
    elif(channel == "2p->ks"):
        (n0, l0, l1) = (2, 1, 0)
    elif final == '2p->kd':
        (n0, l0, l1) = (2, 1, 2)
    elif final == '3d->kp':
        (n0, l0, l1) = (3, 2, 1)
    elif final == '3d->kf':
        (n0, l0, l1) = (3, 2, 3)
    else:
        raise Exception("unsuppoerted channel: "+channel)
        
    hatom = l2.HAtom(1.0)
    driven_term = hatom.length(n0, l0, l1)
    l_op = hatom.h_minus_ene_op(l1, ene)
        
    return (l_op, driven_term)

def get_with_default(dict_obj, key, default):
    if(key in dict_obj):
        return dict_obj[key]
    else:
        return default

# ==== Newton Method ====
def geometric_grad(grad, ab, num_et = None):
    if(num_et == None):
        return geometric_grad_full(grad, ab)
    else:
        return geometric_grad_part(grad, ab, num_et)

def geometric_grad_full(grad, ab):
    a = ab[0]  # first element of geometric sequence
    b = ab[1]  # ratio of geometric sequence
    ns = range(len(grad))
    dxi_da = [b**n for n in ns]        # dxi/da
    dxi_db = [n*a*b**(n-1) for n in ns]  # dxi/db
    da = sum([g*a for (g, a) in zip(grad, dxi_da)])
    db = sum([g*b for (g, b) in zip(grad, dxi_db)])
    return np.array([da, db])

def geometric_grad_part(grad, abxs, num_et):
    ab = [abxs[0], abxs[1]]

    et_index_list = range(num_et)
    et_grad = geometric_grad_full([grad[i] for i in et_index_list], ab)

    ind_index_list = [i for i in range(len(grad)) if i not in et_index_list]
    ind_grad = [grad[i] for i in ind_index_list]

    return np.array(list(et_grad) + list(ind_grad))

def geometric_hess(grad, hess, ab, num_et = None):
    if(num_et == None):
        return geometric_hess_full(grad, hess, ab)
    else:
        return geometric_hess_part(grad, hess, ab, num_et)
    
def geometric_hess_full(grad, hess, ab):
    a = ab[0]
    b = ab[1]
    ns = range(len(grad))
    dx_da    = [b**n for n in ns]                # {dxi/da}_i
    dx_db    = [n*a*b**(n-1) for n in ns]        # {dxi/db}_i
    d2x_da2  = [0 for n in ns]                   # {d2xi/da2}_i
    d2x_dadb = [n*b**(n-1) for n in ns]          # {d2xi/dadb}_i
    d2x_db2  = [a*n*(n-1)*b**(n-2) for n in ns]  # {d2xi/db2}_i
    
    da2 = np.dot(dx_da, np.dot(hess, dx_da)) + np.dot(grad, d2x_da2)
    dadb= np.dot(dx_da, np.dot(hess, dx_db)) + np.dot(grad, d2x_dadb)
    db2 = np.dot(dx_db, np.dot(hess, dx_db)) + np.dot(grad, d2x_db2)

    return np.array([[da2, dadb],
                     [dadb, db2]])


def geometric_hess_part(grad, hess, abxs, num_et):
    ab = [abxs[0], abxs[1]]

    num = len(grad)
    
    et_index_list = range(num_et)
    et_g = grad[0:num_et]
    et_h = hess[0:num_et, 0:num_et]    
    et_hess = geometric_hess_full(et_g, et_h, ab)

    ind_index_list = [i for i in range(len(grad)) if i not in et_index_list]
    ind_hess = hess[num_et:num, num_et:num]

    dxi_da = np.array([ab[1]**n for n in range(num_et)])          # dxi/da
    dxi_db = np.array([n*ab[0]*ab[1]**(n-1) for n in range(num_et)])  # dxi/db

    ind_d2_dxda = np.dot(hess[num_et:num,0:num_et], dxi_da)
    ind_d2_dxdb = np.dot(hess[num_et:num,0:num_et], dxi_db)

    ind_et = np.array([ind_d2_dxda,
                       ind_d2_dxdb])

    return np.r_[np.c_[et_hess,  ind_et],
                 np.c_[ind_et.T, ind_hess]]

def newton(f, z0s, iter_max = 30, eps = 0.0000001, show_lvl = 0):
    """ Newton method for multivariable real/complex function.
    
    Inputs
    ------
    f: lambda [scalar] -> (scalar, [scalar], [[scalar]])
       compute function value, its gradient and its hessians
    z0s: [scalar]
       initial guess
    iter_max: int (30)
       maximum iteration
    eps: double (0.0000001)
       convergence epsilon
    show_lvl: int (0)
       0 => print nothing
       1 => print iteration number i and zi, f(zi), abs_oo(grad)

    Outputs
    -------
    results: bool,     convergence or not
    zs     : [scalar], convergence point
    f(zs)  : scalar,   function value at zs
    grad   : [scalar], gradient at zs
    """

    zs = z0s
    for i in range(iter_max):
        (val, grad, hess)  = f(zs)
        if(show_lvl == 1):
            print i, val, max(abs(grad))

        dz = -la.solve(hess, grad)
        if(max(abs(grad)) < eps and max(abs(dz)) < eps):
            return (True, zs, val, grad)
        zs += dz
    return (False, zs, val, grad)

# ====== Solve only linear program ======
def solve(basis_set, driven_term = None, lop = None, channel = None, energy = None):

    if(channel != None and energy != None):
        (lop, driven_term) = channel_to_op_and_driv(channel, energy)

    if(lop == None and driven_term == None):
        msg = "usage: solve(basis_set, driv, lop) or "
        msg += "solve(basis_set,channel, energy)"
        
        raise Exception()

    def convert(basis):
        if type(basis) == l2.STO:
            return l2.normalized_sto(basis.n, basis.z)
        elif type(basis) == l2.GTO:
            return l2.normalized_gto(basis.n, basis.z)
        else:
            return basis
    
    A00 = np.array([[ l2.cip(a, lop, b) for a in basis_set] for b in basis_set])
    a0  = np.array( [ l2.cip(a, driven_term) for a in basis_set] )
    coefs = la.solve(A00, a0)
    alpha = np.dot(a0, coefs)
    psi = l2.linear_combination(coefs, basis_set)
    return (psi, coefs, alpha)

# ====== Optimize =======
def val_grad_hess(basis_set, driven_term, l_op, opt_index=None):
    """ compute (value,grad,hess) of aA^{-1}a"""

    def make_vec(flag_basis_set, ip):
        return np.array([ip(u, driven_term) if f else 0
                         for (f, u) in flag_basis_set])

    def make_mat(flag_basis_set1, flag_basis_set2, ip):
        return np.array([[ip(a, l_op, b) if fa and fb else 0
                          for (fb, b) in flag_basis_set2]
                         for (fa, a) in flag_basis_set1])

    if opt_index == None:
        opt_index = range(len(basis_set))

    bt = type(basis_set[0])

    us    = [(True, l2.d0_basis(bt, u.n, u.z)) 
             for (u, i) in with_index(basis_set)]
    d_us  = [(i in opt_index, l2.d1_basis(bt, u.n, u.z))
             for (u, i) in with_index(basis_set)]
    dd_us = [(i in opt_index, l2.d2_basis(bt, u.n, u.z))
             for (u, i) in with_index(basis_set)]

    a0 = make_vec(us, l2.cip)
    A00 = make_mat(us, us, l2.cip)

    a1 = make_vec(d_us, l2.cip)
    A10 = make_mat(d_us, us, l2.cip)
    
    a2 = make_vec(dd_us, l2.cip)
    A11 = make_mat(d_us, d_us, l2.cip)
    A20 = make_mat(dd_us, us, l2.cip)

    datas = { 'A00': A00, 'A10':A10, 'A11':A11, 'A20':A20, 
              'a0':a0, 'a1':a1, 'a2':a2 }
    (val, grad, hess) = l_algebra.calc_val_grad_hess_aAm1a(a0, a1, a2, A00, A10, A20, A11)
    return (val, grad[opt_index], hess[opt_index].T[opt_index], datas)


def optimize_simple(init_basis_set, driven_term, l_op, **keyargs):

    res_list = []
    ns = [u.n for u in init_basis_set]
    z0s= [u.z for u in init_basis_set]
    bt = type(init_basis_set[0])

    def one(zs):
        us = [bt(1.0, n, z) for (n, z) in zip(ns, zs)]
        (val, grad, hess, data) = val_grad_hess(us, driven_term, l_op)
        return (val, grad, hess)

    (res_conv, zs_opt, val_opt, grad_opt) = newton(one, z0s, **keyargs)
    basis_opt = [l2.d0_basis(bt, n, z) for (n, z) in zip(ns, zs_opt)]
    return (res_conv, basis_opt, res_list)


def optimize_partial(init_basis_set, driven_term, l_op, 
                     eps=0.00001, max_iter = 30, opt_index=None):
    if(opt_index == None):
        return optimize_partial(init_basis_set, driven_term, l_op, 
                                eps=eps, max_iter = max_iter)
    
    conv_q = convergence_q(eps)
    res_list = []
    basis_set = init_basis_set
    for i in range(max_iter):

        # compute value, gradient, hessian and other
        (v, g, h, data) = val_grad_hess(basis_set, driven_term, l_op, 
                                        opt_index = opt_index)

        # update
        dz_list = la.solve(h, g)
        dz_list_full = list_indexed(zip(opt_index, dz_list), len(basis_set))
        basis_set = [ type(basis_set[0])(1.0, basis.n, basis.z - dz)
                      for (basis, dz) in zip(basis_set, dz_list_full)]
        if conv_q(g, dz_list):
            return (True, basis_set, i+1, res_list)

    return (False, basis_set, max_iter, max_iter, res_list)


def optimize(init_basis_set, driven_term, l_op, 
             eps = 0.00001, 
             iter_callable = None,
             max_iter = 30,
             restriction_index = None,
             opt_index = None):
    """
    restriction_index = [[ 0, 1, 2], [3, 4], [5] ] means that do optimization under restriction with variables
    y0 = x0 = x1 = x2
    y1 = x3 = x4
    y2 = x5

    opt_index = [0, 1] means only optimize of first and second orbital 
    """
    if(opt_index == None):
        opt_index = range(len(init_basis_set))

    conv_q = convergence_q(eps)
    res_list = []
    basis_set = init_basis_set
    for i in range(max_iter):

        # compute value, gradient, hessian and other
        (v, g, h, data) = val_grad_hess(basis_set, driven_term, l_op, 
                                        opt_index = opt_index)

        # modify grad and hessian for restriction if necessary
        if restriction_index == None:
            g_restriction = g
            h_restriction = h
        else:
            g_restriction = [ sum([ g[i] for i in index_list])
                          for index_list in restriction_index]
            h_restriction = [[ sum([ h[i][j] for i in i_list for j in j_list])
                               for i_list in restriction_index]
                             for j_list in restriction_index]

        # update
        dz_list_restriction = la.solve(h_restriction, g_restriction)
        if restriction_index == None:
            i_dz_list = [ (i, dz) for (i, dz) 
                          in zip(range(len(dz_list_restriction)), 
                                 dz_list_restriction)]
        else:
            i_dz_list = [ (i, dz) 
                          for (idx_list, dz) 
                          in zip(restriction_index, dz_list_restriction)
                          for i in idx_list]
        dz_list = [ dz for (i, dz) in sorted(i_dz_list, key = lambda x:x[0])]
        dz_list_full = list_indexed(zip(opt_index, dz_list), len(basis_set))
        basis_set = [ type(basis_set[0])(1.0, basis.n, basis.z - dz)
                      for (basis, dz) in zip(basis_set, dz_list_full)]
        if iter_callable != None:
            res_list.append(iter_callable(basis_set))
        if conv_q(g_restriction, dz_list):
            return (True, basis_set, res_list)

    return (False, basis_set, res_list)

def optimize_hydrogen_pi(basis_set, channel, energy, 
                         eps = 0.00001, 
                         iter_callable = None, 
                         max_iter = 50):
    (l_op, driven_term) = channel_to_op_and_driv(channel, energy)
    return optimize(basis_set, driven_term, l_op, 
                    eps, 
                    iter_callable,
                    max_iter)


def optimize_even_temp(bt, et_num, et_n, et_z0, et_r0, nz0_list, driv, lop, **keyargs):
    """
    Inputs
    ------
    et_num : int
        number of the even tempered basis
    et_n: int
        princinple number of the even tempered basis
    et_z0: complex
        initial guess for first complex orbital exponents of the even-tempered basis
    et_r: complex
        initial guess for ratio of the even tempered basis
    nz_list: [(int, complex)]
        principle numbers and initial guess orbital exponents of basis to be optimize
    """

    res_list = []
    def one(vars):
        a = vars[0]
        b = vars[1]
        et_ns = range(et_num)
        zs = [a*b**n for n in et_ns] + list(vars[2:])
        ns = [et_n for n in et_ns] + [n for (n, z) in nz0_list]
        us = [bt(1.0, n, z) for (n, z) in zip(ns, zs)]
        (v0, g0, h0, data) = val_grad_hess(us, driv, lop)
        grad = geometric_grad(g0, [a, b], et_num)
        hess = geometric_hess(g0, h0, [a, b], et_num)
        return (v0, grad, hess)

    guess = [et_z0, et_r0] + [z for (n, z) in nz0_list]
    (convq, zs_o, val_o, grad_o) = newton(one, guess, **keyargs)
    a_o = zs_o[0]
    b_o = zs_o[1]
    zs_o= zs_o[2:]
    n_o = [n for (n, z) in nz0_list]
    nz_o = [(et_n, a_o*b_o**k) for k in range(et_num)] +\
           [(n, z) for (n, z) in zip(n_o, zs_o)]
    basis_o = [l2.d0_basis(bt, n, z) for (n, z) in nz_o]
    return (convq, basis_o, res_list)

# ====== Optimize for energies ========
def optimize_for_energies_2pks(basis_set0, energy_list,  **keyargs):
    b_type = type(basis_set0[0])
    if b_type !=l2.STO and b_type != l2.GTO:
        print 'unsupported basis'
        sys.exit()
    driven_term = l2.h_like_atom('2p').dipole_init_velocity(0)
    def one_iter(basis_set, energy):
        if b_type == l2.STO:
            l_op = l2.h_like_atom('1s').h_minus_energy_sto(energy)
        else:
            l_op = l2.h_like_atom('1s').h_minus_energy_gto(energy)
        return optimize(basis_set, driven_term, l_op, **keyargs)

    ene_zeta_list = []
    basis_set = basis_set0
    for energy in energy_list:
        (conv, basis_set, dummy) = one_iter(basis_set, energy)
        if not conv:
            return (False, ene_zeta_list)
        zeta_list = [ basis.z for basis in basis_set]
        ene_zeta_list.append( (energy, zeta_list) )
    return (True, ene_zeta_list)


# ===== Phase Shift ======
# see 2015/9/7_cbf_wkb.ipynb
def coulomb_phase_shift(k, L):
    return cmath.phase(gamma(1.0 + L - 1.0j/k))

def wkb_local_momentum(k, L, r):
    c = 1.0 * L * (L + 1)
    return np.sqrt(k**2 + 2.0/r - c/(r**2))

def wkb_zeta(k, L, r):
    p = wkb_local_momentum(k, L, r)
    c = 1.0*L*(L+1)
    term1 = 1.0/(4.0*p**3)*(3.0*c/(r**4)-2.0/(r**3))
    term2 = 5.0/(8.0*p**5)*(1.0/(r*r)-c/(r**3))**2
    return p + term1 + term2

def wkb_zeta_prime(k, L, r):
    L = 1.0 * L
    return -(3.0*L*(L + 1)/r**5 - 3.0/(2*r**4))/(-L*(L + 1)/r**2 + k**2 + 2/r)**(3.0/2.0) + 5.0*(-12*L*(L + 1)/r**4 + 8/r**3)*(2*L*(L + 1)/r**3 - 2/r**2)/(32*(-L*(L + 1)/r**2 + k**2 + 2/r)**(5.0/2.0)) - (-3*L*(L + 1)/(4*r**4) + 1/(2*r**3))*(-3*L*(L + 1)/r**3 + 3/r**2)/(-L*(L + 1)/r**2 + k**2 + 2/r)**(5.0/2.0) + 5*(-5*L*(L + 1)/r**3 + 5/r**2)*(2*L*(L + 1)/r**3 - 2/r**2)**2/(32*(-L*(L + 1)/r**2 + k**2 + 2/r)**(7.0/2.0)) + (L*(L + 1)/r**3 - 1/r**2)/sqrt(-L*(L + 1)/r**2 + k**2 + 2/r)

def wkb_phi(k, L, r):
    p = wkb_local_momentum(k, L, r)
    if L == 0:
        phi0 = p*r + 1.0/k*np.log(2+2*k*k*r+2*p*k*r)
        term1 = (9+6*k**2*r+2*k**4*r**2)/(24*p**3*r**2)
        return phi0 + term1
    else:
        c = 1.0 * L * (L + 1)
        phi0 = p*r - np.sqrt(c)*np.arctan((r-c)/(np.sqrt(c)*r*p)) + 1.0/k*np.log(2+2*k*k*r+2*p*k*r)
        term1 = - (c + 2*c**2*k**2 + 3*c*k**4*r**2 + 3*r + 6*k**2*r**2 + k**4*r**3) / (24.0*r**3*p**3*(1+c*k*k))
        term2 = - 1.0/(8*np.sqrt(c)) * np.arctan((r-c)/(np.sqrt(c)*p*r))
        return phi0 + term1 + term2

def wkb_delta_and_A(k, L, r, psi, d_psi):
    zeta = wkb_zeta(k, L, r)
    zeta_prime = wkb_zeta_prime(k, L, r)
    # log_d_psi = d_psi/psi
    # arg_arctan = zeta/(log_d_psi+zeta_prime/(2.0*zeta))
    arg_arctan = zeta*psi/(d_psi + psi*zeta_prime/(2.0*zeta))
    delta = np.arctan(arg_arctan)
    A = np.sqrt(zeta) * psi / np.sin(delta)
    return (delta, A)

def wkb_func(k, L, r_match, psi, d_psi):
    (delta, A) = wkb_delta_and_A(k, L, r_match, psi, d_psi)
    phi0 = wkb_phi(k, L, r_match)
    def __func__(r):
        phi = wkb_phi(k, L, r) - phi0
        zeta = wkb_zeta(k, L, r)
        return A * zeta**(-0.5) * np.sin(phi+delta)
    return __func__

def wkb_phase_shift(k, L, r, psi, d_psi, r_asym):
    (delta, A) = wkb_delta_and_A(k, L, r, psi, d_psi)
    phi = wkb_phi(k, L, r_asym) - wkb_phi(k, L, r)
    phase_shift = phi + delta - k*r_asym - (1.0/k)*np.log(2.0*k*r_asym) + L*np.pi/2.0

    if A > 0.0:
        return phase_shift % (2*np.pi)
    else:
        return (phase_shift + np.pi) % (2*np.pi)
    
def phase_shift_by_coulomb_func(psi_r, r, r_match, k):
    """
    psi_r: caluculated wave function written in sympy symbol r
    r:     sympy symbol used in psi_r
    r_match: matching point 
    k:       wave number
    """
    
