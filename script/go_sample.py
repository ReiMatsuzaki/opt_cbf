from cmath import rect
import time
sys.path.append("~/src/git/opt_cbf/script")
from opt_cbf import *

def create_z_list():
    x_list = [ 0.001 * 2.0 ** n for n  in range(10)]
    return [complex(x,-y) for x in x_list for y in x_list]

def create_r_list():
    r_abs_list = [ 0.1 * n for n in [7,8,9,10,11,12,13]]
    r_arg_list = [ 3.1415 * n / 180.0 for n in range(1,20,2)]
    return [rect(r, -phi) for r in r_abs_list for phi in r_arg_list]


f = open("2et_g.in", "r")
str_base =  f.read()
f.close()
z_list = create_z_list()[0:10]
r_list = create_r_list()[0:10]

zz_array = search(str_base, z_list, r_list)
t0 = time.clock()
zz_uniq  = take_uniq(zz_array, near)
t1 = time.clock()
print "uniq_time:" + str(t1-t0)

for zz in zz_uniq:
    print zz


