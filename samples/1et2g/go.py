from cmath import rect
import time
import sys
import os

sys.path.append(os.environ["OPT_CBF"] + "/script")
from opt_cbf import *

opt_cbf  = os.environ["OPT_CBF"] + "/opt_cbf tmp.in tmp.out"
print opt_cbf

def create_z_list():
    x_list = [ 0.001 * 2.0 ** n for n  in range(10)]
    return [complex(x,-y) for x in x_list for y in x_list]

def create_r_list():
    r_abs_list = [ 0.1 * n for n in [7,8,9,10,11,12,13]]
    r_arg_list = [ 3.1415 * n / 180.0 for n in range(1,20,2)]
    return [rect(r, -phi) for r in r_abs_list for phi in r_arg_list]

f = open("base.in", "r")
str_base =  f.read()
f.close()
z_list = create_z_list()[0:10]
r_list = create_r_list()[0:10]

t0_calc = time.clock()
zz_array = []
for z in z_list:
    for r in r_list:
        create_in_et(z, r, str_base, 'tmp.in')
        os.system(opt_cbf)
        out_file = open('tmp.out')
        strs_out = out_file.readlines()
        out_file.close()
        if(ok_conv(strs_out) and 
           ok_coef(strs_out, 0.001) and 
           ok_et_basis(strs_out, 0.001)):
            zz_array.append( data_et(strs_out)[0] )
t1_calc = time.clock()

t0_uniq = time.clock()
zz_uniq  = take_uniq(zz_array, near_et(0.001))
t1_uniq = time.clock()

print "number:" + str(len(z_list) * len(r_list))
print "calc_time:" + str(t1_calc-t0_calc)
print "uniq_time:" + str(t1_uniq-t0_uniq)

for zz in zz_uniq:
    print zz


