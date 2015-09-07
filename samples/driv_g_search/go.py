import sys
import os
sys.path.append("/Users/rei/src/git/var_param/")
import var_param as vp
opt_cbf_root = "/Users/rei/src/git/opt_cbf"
sys.path.append(opt_cbf_root + "/script/")
import opt_cbf as oc
import time

opt_cbf_in_file = "opt_cbf.in"
opt_cbf_out_file = "opt_cbf.out"

def convergence():
    strs_out = open(opt_cbf_out_file).readlines()
    res = oc.ok_conv(strs_out) and oc.ok_coef(strs_out, 0.001) and oc.ok_opt_basis(strs_out, 0.001)
    return res

cmd = "{0}/opt_cbf {1} {2}".format(opt_cbf_root,
                                   opt_cbf_in_file, 
                                   opt_cbf_out_file)

file_list = [ opt_cbf_in_file ]
xs = [ 0.001 * 2.0 **n for n in range(10) ]
print xs

t0 = time.clock()
vp.run(
    root = os.path.dirname(__file__),
    files = file_list,
    params = vp.params.combinations(["x0r","x0i","x1r","x1i","x2r","x2i"], xs),
    commands = vp.command(cmd),
    out_dir = vp.out_dir.flat_from("results"),
    save_if = convergence
)
t1 = time.clock()
print "timer: ", t1-t0

    
