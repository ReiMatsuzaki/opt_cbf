import sys
import os
import os.path
import numpy as np
import matplotlib.pyplot as plt

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

    print idx

    # directory
    dir = str(idx)
    if not os.path.exists(dir):
        os.mkdir(dir) 
    os.chdir(dir)

    # extract input data
    d = [l.strip() for l in line[1:-2].split(",")]
    z=complex(d[2])
    r=complex(d[3])

    # create input file and calculation
    create_in_et(z, r, str_base, 'opt_cbf.in')
    os.system(opt_cbf_path +  '/opt_cbf opt_cbf.in opt_cbf.out')

    # load wave functions data and plot it as 
    psi_data = np.loadtxt('psi.out', delimiter=',')
    r  = psi_data[:,0] 
    re = psi_data[:,1]
    im = psi_data[:,2]
    plt.plot(r, re)
    plt.plot(r, im)
    plt.savefig('psi.eps')
    plt.clf()   # clear figure

    # end
    os.chdir('..')
    idx += 1

os.chdir('..')
search_out.close()

