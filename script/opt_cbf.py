import time
import os
import re

# 
# ===================================
#  Utilities
def take_uniq_one(xs_ys, eq) :
    (xs, ys) = xs_ys
    if(ys == None or len(ys) == 0):
        return xs
    elif(len(ys) == 1):
        return xs + [ys[0]]
    y0   = ys[0]
    rest = [y for y in ys if not(eq(y, y0))]
    return take_uniq_one((xs + [y0], rest), eq)

def take_uniq(xs, eq): 
    return take_uniq_one(([], xs), eq)

def take_and(xs):
    return reduce(lambda a,b : a and b, xs)


#
# ==================================
#  cast
#
# string->complex
def cast_to_complex(line):
    m = re.search('\((.*), *(.*)\)', line)
    if(m == None):
        print "failed to convert in complex_cpp2py"
    z = complex(float(m.group(1)), float(m.group(2)))
    return z

# string->bool
def cast_to_bool(str):
    str2 = str.strip().lower()
    if str2 == "true":
        return True
    elif str2 == "false":
        return False
    else:
        return None

# [string->?]->(string->[?])
def cast_array_to_array(s0, fs):
    ss = s0.split(" ")
    return [f(s) for (f,s) in zip(fs, ss)]

def cast_array(fs):
    return lambda s: cast_array_to_array(s, fs)

# string -> (complex,complex)
def cast_iicc_to_cc(str):
    return tuple(map(cast_to_complex, str.strip().split(" ")[2:4]))


# 
# ==================================
# [string],string,(string->?)->[?]
# (str_keys_values) find values whose key is /key/
def values_for_key(str_keys_values, key, cast=None):

    vals = [line.split(":")[1].strip()
            for line 
            in str_keys_values
            if line.count(key)]
    if(cast == None):
        return vals
    else:
        if(type(cast) == list):
            f = cast_array(cast)
        else:
            f = cast

        return [f(v) for v in vals]

#
# ==================================
#  Check validity of resutls
# from [string] created by readlines method, return is 
# valid output of opt_cbf program
#
def ok_conv(strs_out):
    return values_for_key(strs_out,
                          "convergence",
                          cast_to_bool)[0]

def ok_coef(strs_out, eps):
    cs = values_for_key(strs_out, "coef", cast_to_complex)
    c0 = min([abs(c) for c in cs])
    return c0 > eps

def ok_et(strs_out, eps):
    zr_list = values_for_key(strs_out, "opt_et_basis",
                                cast_iicc_to_cc)
    ok_zeta = [z.real > 0 and z.imag < 0 and r.imag < 0
               for (z,r) in zr_list]
    ok_ratio = [ abs(r-1) > eps for (z,r) in zr_list]
    return take_and(ok_zeta + ok_ratio)

def ok_opt_basis(strs_out, eps):
    zs = zetas_opt_cbf(strs_out);
    
    ok_dim = [z.real > 0 and z.imag < 0 for z in zs]

    num = len(zs)
    dists = [abs(zs[i] - zs[j]) for i in range(num)
             for j in range(i)]
    ok_dist = [min(dists) > eps]

    return take_and(ok_dim + ok_dist)

#
# ==================================
# ET-Basis optimization
#
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


#
# ==================================
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
    
def uniq_opt_cbf(zs_list, eps):
    return take_uniq(zs_list, near_opt_cbf(eps))
