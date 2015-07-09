from cmath import rect
import time
import sys
sys.path.append("/Users/rei/src/git/opt_cbf/script")
from opt_cbf import *

def create_z_list():
    x_list = [ 0.001 * 2.0 ** n for n  in range(10)]
    return [complex(x,-y) for x in x_list for y in x_list]

f = open("base.in", "r")
str_base =  f.read()
f.close()
zs = create_z_list()[0:10]

t0 = time.clock()
zs_array = []
num = len(zs)
for i in range(num):
    for j in range(i):
        z0 = zs[i]
        z1 = zs[j]
        create_in_opt_cbf([z0,z1], str_base, 'tmp.in')
        os.system('${HOME}/src/git/opt_cbf/opt_cbf tmp.in tmp.out')
        out_file = open('tmp.out')
        out_kv   = out_file.readlines()
        out_file.close()
        if(ok_conv(out_kv) and ok_coef(out_kv, 0.001) and
           ok_opt_basis(out_kv, 0.001)):
            zetas = zetas_opt_cbf(out_kv)
            zs_array.append(zetas)

t1 = time.clock()
print "calc_time:" + str(t1-t0)

zs_array_uniq = uniq_opt_cbf(zs_array, 0.001)

for zs in zs_array_uniq:
    print zs


