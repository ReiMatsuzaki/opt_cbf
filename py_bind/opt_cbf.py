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
        
def channel_to_op_and_driv(channel, ene, basis_type):

    [init, final] = channel.split("->")    
    if final == 'ks':
        l1 = 0
    elif final == 'kp':
        l1 = 1
    elif final == 'kd':
        l1 = 2
    elif final == 'kf':
        l1 = 3
        
    h_init = l2.h_like_atom(init)    
    driven_term = h_init.dipole_init_length(l1)
    h_final = l2.h_like_atom(str(l1 + 1) + final[1])
    if basis_type == l2.STO:
        l_op = h_final.h_minus_energy_sto(ene)
    elif basis_type == l2.GTO:
        l_op = h_final.h_minus_energy_gto(ene)
    else:
        print 'opt_cbf.py::channel_to_op_and_driv'
        print 'unsupported basis_type: ' + basis_type
        
    return (l_op, driven_term)


# ====== Solve only linear program ======
def solve(basis_set, driven_term, lop):

    def convert(basis):
        if type(basis) == l2.STO:
            return l2.normalized_sto(basis.n, basis.z)
        elif type(basis) == l2.GTO:
            return l2.normalized_gto(basis.n, basis.z)
        else:
            return basis
    
    A00 = np.array([[ l2.cip(a, l_op(b)) for a in us] for b in basis_set])
    a0  = np.array( [ l2.cip(a, driven_term) for a in basis_set] )
    coefs = la.solve(A00, a0)
    alpha = np.dot(a0, coefs)
    psi = linear_combination(coefs, basis_set)
    return (psi, coefs, alpha)

# ====== Optimize =======
def val_grad_hess(basis_set, driven_term, l_op, opt_index=None):
    """ compute (value,grad,hess) of aA^{-1}a"""

    ip_mat = get_sym_ip(type(basis_set[0]), type(basis_set[0]))
    ip_vec = get_sym_ip(type(basis_set[0]), type(driven_term.prim_i(0)))

    if type(basis_set[0]) == l2.STO:
        normalized    = l2.normalized_sto
        d_normalized  = l2.d_normalized_sto
        dd_normalized = l2.dd_normalized_sto
    elif type(basis_set[0]) == l2.GTO:
        normalized    = l2.normalized_gto
        d_normalized  = l2.d_normalized_gto
        dd_normalized = l2.dd_normalized_gto
    else:
        print 'opt_cbf.py/val_grad_hess'
        print 'unsupported type: ' + type(basis_type[0])
        sys.exit(0)

    def make_vec(flag_basis_set):
        return np.array([ip_vec(u, driven_term) if f else 0
                         for (f, u) in flag_basis_set])

    def make_mat(flag_basis_set1, flag_basis_set2):
        return np.array([[ip_mat(a, l_op(b)) if fa and fb else 0
                          for (fb, b) in flag_basis_set2]
                         for (fa, a) in flag_basis_set1])

    if opt_index == None:
        opt_index = range(len(basis_set))

    us    = [(True, normalized(u.n, u.z)) 
             for (u, i) in with_index(basis_set)]
    d_us  = [(i in opt_index, d_normalized(u.n, u.z))
             for (u, i) in with_index(basis_set)]
    dd_us = [(i in opt_index, dd_normalized(u.n, u.z))
             for (u, i) in with_index(basis_set)]

    a0 = make_vec(us)
    A00 = make_mat(us, us)

    a1 = make_vec(d_us)
    A10 = make_mat(d_us, us)
    
    a2 = make_vec(dd_us)
    A11 = make_mat(d_us, d_us)
    A20 = make_mat(dd_us, us)

    datas = { 'A00': A00, 'A10':A10, 'A11':A11, 'A20':A20, 
              'a0':a0, 'a1':a1, 'a2':a2 }
    (val, grad, hess) = l_algebra.calc_val_grad_hess_aAm1a(a0, a1, a2, A00, A10, A20, A11)
    return (val, grad[opt_index], hess[opt_index].T[opt_index], datas)


def optimize_simple(init_basis_set, driven_term, l_op, 
                    eps = 0.00001, max_iter = 30):
    conv_q = convergence_q(eps)
    res_list = []
    basis_set = init_basis_set
    for i in range(max_iter):

        # compute value, gradient, hessian and other
        (v, g, h, data) = val_grad_hess(basis_set, driven_term, l_op)

        # update
        dz_list = la.solve(h, g)
        basis_set = [ type(basis_set[0])(1.0, basis.n, basis.z - dz)
                      for (basis, dz) in zip(basis_set, dz_list)]
        if conv_q(g, dz_list):
            return (True, basis_set, res_list)

    return (False, basis_set, res_list)

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

    opt_index = [0, 1] : only optimize of first and second orbital 
    """

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

def optimize_hydrogen_pi(basis_set, channel, ene, 
                         eps = 0.00001, 
                         iter_callable = None, 
                         max_iter = 50):
    (l_op, driven_term) = channel_to_op_and_driv(channel, ene, type(basis_set[0]))
    return optimize(basis_set, driven_term, l_op, 
                    eps, 
                    iter_callable,
                    max_iter)


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
    
