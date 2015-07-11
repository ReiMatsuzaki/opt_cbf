import sys
import os
import os.path

opt_cbf_path = os.environ["OPT_CBF"] 
script_path = opt_cbf_path + "/script"
print script_path
sys.path.append(script_path)
from opt_cbf import *

# output file for search 
search_out = open('search.out')

# skip unnecessary data line
for i in range(4):
    search_out.readline()

# copy base file as string
f = open('base_write.in')
str_base = f.read()
f.close()

# read necessary datas as [string]
lines = search_out.readlines()

base_dir = 'each_calc'
if not os.path.exists(base_dir):
    os.mkdir(base_dir)
os.chdir(base_dir)

idx = 0
for line in lines:
    d = [l.strip() for l in line[1:-2].split(",")]
    z=complex(d[2])
    r=complex(d[3])
    dir = str(idx)
    if not os.path.exists(dir):
        os.mkdir(dir) 
    os.chdir(dir)
    create_in_et(z, r, str_base, 'opt_cbf.in')
    print idx
    os.system(opt_cbf_path +  '/opt_cbf opt_cbf.in opt_cbf.out')
    os.chdir('..')
    idx += 1

os.chdir('..')
search_out.close()

