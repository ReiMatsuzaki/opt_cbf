import time
import os
import re

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

def create_in_file(z, r, str_base_file, file_name) :
    l = str_base_file
    l = l.replace("__xr__", str(z.real))
    l = l.replace("__xi__", str(z.imag))
    l = l.replace("__rr__", str(r.real))
    l = l.replace("__ri__", str(r.imag))
    
    f = open(file_name, "w")    
    f.write(l)
    f.close()

# from string array (str_keys_values) find values whose key is /key/
def values_for_key(str_keys_values, key, cast_func=None):
    vals = [line.split(":")[1].strip()
            for line 
            in str_keys_values
            if line.count(key)]
    if(cast_func == None):
        return vals
    else:
        return [cast_func(v) for v in vals]

def take_and(xs):
    return reduce(lambda a,b : a and b, xs)

def is_convergenceQ(strs_kv):

    # check for convergence
    lines = strs_kv
    ok_convergence = values_for_key(lines, "convergence", cast_to_bool)[0]

    # check for non zero coefficient
    coef_list = values_for_key(lines, "coef", cast_to_complex)
    min_coef_abs = min(map(abs, coef_list))
    ok_coef = min_coef_abs > 0.001

    # check for 3rd region in complex
    z0_r0_list = values_for_key(strs_kv, "opt_et_basis", cast_iicc_to_cc)
    
    ok_3rd = take_and(
        [(z0.real > 0 and z0.imag < 0 and r0.imag < 0)
         for (z0,r0) in z0_r0_list])

    # check ratio is not near one
    ok_ratio = take_and([ abs(r0-1) > 0.001 for (z0,r0) in z0_r0_list])

    return take_and([ok_convergence, ok_coef, ok_3rd, ok_ratio])
  
def cast_to_complex(line):
    m = re.search('\((.*), *(.*)\)', line)
    if(m == None):
        print "failed to convert in complex_cpp2py"
    z = complex(float(m.group(1)), float(m.group(2)))
    return z

def cast_to_bool(str):
    str2 = str.strip().lower()
    if str2 == "true":
        return True
    elif str2 == "false":
        return False
    else:
        return None

# string of (int,int,complex,complex)->(complex,complex)
def cast_iicc_to_cc(str):
    return tuple(map(cast_to_complex, str.strip().split(" ")[2:4]))

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
                zz_array.append(values_for_key(out_kv, "opt_et_basis", cast_iicc_to_cc))
    t1 = time.clock()
    print "Time: " + str(t1-t0)
    return zz_array

def near (a, b, eps):
    f = lambda x,y: abs(x,y)<eps
    take_and([ f(a0,b0) and f(ar,br) for ((a0,ar),(b0,br)) in zip(a,b)])
