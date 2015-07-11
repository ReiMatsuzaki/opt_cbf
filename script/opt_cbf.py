import time
import os
import re

# ====== Utilities ==================

# take_uniq extract non equal elements from list xs
# the equality is calculated from eq.
# take_uniq_one is one process for take_uniq
# ([a], [a]) -> (a,a->bool)
def take_uniq_one(xs_ys, eq) :
    (xs, ys) = xs_ys
    if(ys == None or len(ys) == 0):
        return xs
    elif(len(ys) == 1):
        return xs + [ys[0]]
    y0   = ys[0]
    rest = [y for y in ys if not(eq(y, y0))]
    return take_uniq_one((xs + [y0], rest), eq)

# [a],(a,a->bool)->[a]
def take_uniq(xs, eq): 
    return take_uniq_one(([], xs), eq)

# [bool]->bool
def take_and(xs):
    return reduce(lambda a,b : a and b, xs)


# ======= Cast =======================

# string->complex
def cast_to_complex(line):
    m = re.search('\((.*), *(.*)\)', line)
    if(m == None):
        print "failed to convert in complex_cpp2py"
    z = complex(float(m.group(1)), float(m.group(2)))
    return z

# string->bool
def cast_to_bool(s):
    s1 = s.strip().lower()
    if s1 == "true":
        return True
    elif s1 == "false":
        return False
    else:
        return None

# string,[string->?]->[?]
def cast_array_to_array(s0, fs):
    ss = s0.split(" ")
    return [f(s) for (f,s) in zip(fs, ss)]

# [string->?] -> (string->[?])
def cast_array(fs):
    return lambda s: cast_array_to_array(s, fs)

# [string],string,(string->?)->[?]
# from kvs find values whose key is k
def values_for_key(kvs, k, cast=None):

    vals = [line.split(":")[1].strip()
            for line 
            in kvs
            if line.count(k)]
    if(cast == None):
        return vals
    else:
        if(type(cast) == list):
            f = cast_array(cast)
        else:
            f = cast

        return [f(v) for v in vals]


# ======== Checker ======================
#  Check validity of resutls
# from [string] created by readlines method, return is 
# valid output of opt_cbf program
# [string] -> bool
def ok_conv(ss):
    return values_for_key(ss,
                          "convergence",
                          cast_to_bool)[0]

def ok_coef(ss, eps):
    cs = values_for_key(ss, "coef", cast_to_complex)
    c0 = min([abs(c) for c in cs])
    return c0 > eps

def ok_et_basis(ss, eps):
    data_list = data_et(ss)

    zs = [ z * r**k for (n,num,z,r) in data_list for k in range(num) ]
    ok_zeta = [ z.imag < 0 and z.real > 0 for z in zs]

    ok_ratio = [ r.imag < 0 for (n,num,z,r) in data_list]
                 
    ok_ratio2 = [ abs(r-1) > eps for (n,num,z,r) in data_list]

    return take_and(ok_zeta + ok_ratio + ok_ratio2)

def ok_opt_basis(ss, eps):
    zs = zetas_opt_cbf(ss);
    
    ok_dim = [z.real > 0 and z.imag < 0 for z in zs]

    num = len(zs)
    dists = [abs(zs[i] - zs[j]) for i in range(num)
             for j in range(i)]
    ok_dist = [min(dists) > eps]

    return take_and(ok_dim + ok_dist)


# ========== ET-Basis ==============
# to be rmeoved
def search(str_base, z_list, r_list):
    print "NumSearch: " + str(len(z_list) * len(r_list) )
    t0 = time.clock()
    zz_array = []
    for z in z_list:
        for r in r_list:
            create_in_file(z, r, str_base, "tmp.in")
            os.system('${HOME}/bin/opt_cbf tmp.in tmp.out')
            out_file = open('tmp.out')
            out_kv = out_file.readlines()
            out_file.close()
            if(is_convergenceQ(out_kv)):
                iicc = values_for_key(out_kv, "opt_et_basis", [int,int,cast_to_complex, cast_to_complex])
                cc = (iicc[2], iicc[3])
                zz_array.append(cc)
    t1 = time.clock()
    print "Time: " + str(t1-t0)
    return zz_array

# to be removed
def near (a, b, eps):
    f = lambda x,y: abs(x,y)<eps
    take_and([ f(a0,b0) and f(ar,br) for ((a0,ar),(b0,br)) in zip(a,b)])

def create_in_et(z, r, str_base_file, file_name) :
    l = str_base_file
    l = l.replace("__xr__", str(z.real))
    l = l.replace("__xi__", str(z.imag))
    l = l.replace("__rr__", str(r.real))
    l = l.replace("__ri__", str(r.imag))
    
    f = open(file_name, "w")    
    f.write(l)
    f.close()

def create_in_multi_et(zs, rs, str_base_file, file_name):
    l = str_base_file
    num = len(zs)
    for i in range(num):
        l = l.replace("__x%sr__"%i, str(zs[i].real))
        l = l.replace("__x%si__"%i, str(zs[i].imag))
        l = l.replace("__r%sr__"%i, str(rs[i].real))
        l = l.replace("__r%si__"%i, str(rs[i].imag))

    f = open(file_name, "w")
    f.write(l)
    f.close()
        

def near_et(eps):

    def __near__(x,y):
        (a,b,z0,r0) = x
        (c,d,z1,r1) = y
        return abs(z0-z1) < eps and abs(r0-r1) < eps
    return __near__

# [string]->[(int,int,complex, complex)]
def data_et(ss):
    caster = [int,int,cast_to_complex, cast_to_complex]
    iicc_list = values_for_key(ss, "opt_et_basis", caster)
    return iicc_list


# ========== Opt-Basis =============
# Independenet optimization
#
def create_in_opt_cbf(zs, str_base_file, file_name):

    l = str_base_file
    idx = 0
    for z in zs:
        k_r = "__z" + str(idx) + "r__"
        l=l.replace("__z" + str(idx)+ "r__", str(z.real))
        l=l.replace("__z" + str(idx)+ "i__", str(z.imag))
        idx = idx + 1

    f = open(file_name, "w")
    f.write(l)
    f.close()

# ss_kv :: keys_values list
def zetas_opt_cbf(ss_kv):
    cast = cast_array((int,cast_to_complex))
    ics = values_for_key(ss_kv, "opt_basis", cast)
    zetas = [ c for (i,c) in ics]    
    return zetas

# double -> ([complex], [complex] -> bool)
def near_opt_cbf(eps):

    def __near(xs, ys):
        f = lambda x: (x.real, x.imag)
        xs0 = sorted(xs, key=f)
        ys0 = sorted(ys, key=f)
        dmin = min([ abs(x-y) for (x,y) in zip(xs0, ys0)])
        return (dmin < eps)
    return __near
    
# [complex], double -> [complex]
def uniq_opt_cbf(zs_list, eps):
    return take_uniq(zs_list, near_opt_cbf(eps))
