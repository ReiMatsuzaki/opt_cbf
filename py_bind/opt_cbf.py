import sys
import copy
import numpy as np
import scipy.linalg as la
import l_algebra
import l2func as l2

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
def solve(basis_set, channel = None, energy = None, 
                 driven_term = None, l_op = None):
    if channel != None and energy != None:
        if (type(basis_set[0]) == l2.STO or type(basis_set[0]) == l2.STOs):
            b_type = l2.STO
        elif (type(basis_set[0]) == l2.GTO or type(basis_set[0]) == l2.GTOs):
            b_type = l2.GTO
        else:
            print 'opt_cbf.py/solve'
            print 'unsupported type: ' + type(basis_type[0])
        (l_op, driven_term) = channel_to_op_and_driv(channel, energy, b_type)

    if l_op == None or driven_term == None:
        print "use this function such as..."
        print "solve_by_cbf(basis_set, channel = '1s->kp', energy = 0.5)"
        print "solve_by_cbf(basis_set, driven_term = d_term, l_op = l_op)"
        sys.exit()

    def convert(basis):
        if type(basis) == l2.STO:
            return l2.normalized_sto(basis.n, basis.z)
        elif type(basis) == l2.GTO:
            return l2.normalized_gto(basis.n, basis.z)
        else:
            return basis

    us = [ convert(basis) for basis in basis_set]
    t_basis = type(us[0].prim_i(0))    
    t_driv  = type(driven_term.prim_i(0))
    ip_mat = get_sym_ip(t_basis, t_basis)
    ip_vec = get_sym_ip(t_basis, t_driv)
    
    A00 = np.array([[ ip_mat(a, l_op(b)) for a in us] for b in us])
    a0  = np.array( [ ip_vec(a, driven_term) for a in us] )
    coefs = la.solve(A00, a0)
    alpha = np.dot(a0, coefs)

    cum_psi = type(us[0])() # call constructor of same class 
    for (c, u) in zip(coefs, us):
        cum_psi.add(l2.scalar_prod(c, u))
    return (cum_psi, coefs, alpha)


# ====== Optimize =======
def val_grad_hess(basis_set, driven_term, l_op):
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
        

    def make_basis(func, basis):
        return func(basis.n, basis.z)

    us    = [ make_basis(normalized, basis) for basis in basis_set]
#    print [ type(u) for u in us]
    a0    = np.array( [ ip_vec(a, driven_term) for a in us])
    A00   = np.array( [ [ ip_mat(a, l_op(b)) 
                          for a in us] for b in us])

    d_us  = [ make_basis(d_normalized, basis)  
              for basis in basis_set]
    a1    = np.array( [ ip_vec(a, driven_term) for a in d_us])
    A10   = np.array( [ [ ip_mat(a, l_op(b)) 
                          for a in us] for b in d_us])
    
    dd_us = [ make_basis(dd_normalized, basis) 
              for basis in basis_set]
    a2    = np.array( [ ip_vec(a, driven_term) for a in dd_us])
    A11   = np.array( [ [ ip_mat(a, l_op(b)) 
                          for a in d_us] for b in d_us])
    A20   = np.array( [ [ ip_mat(a, l_op(b))
                          for a in us] for b in dd_us])
    
    datas = { 'A00': A00, 'A10':A10, 'A11':A11, 'A20':A20, 
              'a0':a0, 'a1':a1, 'a2':a2 }
    (val, grad, hess) = l_algebra.calc_val_grad_hess_aAm1a(a0, a1, a2, A00, A10, A20, A11)
    return (val, grad, hess, datas)

def optimize(basis_set, driven_term, l_op, 
             eps = 0.00001, 
             iter_callable = None,
             max_iter = 30,
             restriction_index = None):
    """
    restriction_index = [[ 0, 1, 2], [3, 4], [5] ] means that do optimization under restriction with variables
    y0 = x0 = x1 = x2
    y1 = x3 = x4
    y2 = x5
    """

    conv_q = convergence_q(eps)
    res_list = []
    for i in range(max_iter):

        # compute value, gradient, hessian and other
        (v, g, h, data) = val_grad_hess(basis_set, driven_term, l_op)

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
                          for (idx_list, dz) in zip(restriction_index, dz_list_restriction)
                          for i in idx_list]
        dz_list = [ dz for (i, dz) in sorted(i_dz_list, key = lambda x:x[0])]
        basis_set = [ type(basis_set[0])(1.0, basis.n, basis.z - dz)
                       for (basis, dz) in zip(basis_set, dz_list)]
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
def phase_shift_by_coulomb_func(psi_r, r, r_match, k):
    """
    psi_r: caluculated wave function written in sympy symbol r
    r:     sympy symbol used in psi_r
    r_match: matching point 
    k:       wave number
    """
    
