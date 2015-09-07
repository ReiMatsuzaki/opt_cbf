import sys
import os
opt_cbf_root = "/Users/rei/src/git/opt_cbf"
sys.path.append(opt_cbf_root + "/script")
import opt_cbf as oc
sys.path.append("/Users/rei/src/git/var_param/")
import var_param as vp

res_dir = "./results/"

def read_zeta_from_dir(d):
    f = res_dir + d + "/opt_cbf.out"
    datas = open(f).readlines()
    return oc.zetas_opt_cbf(datas)

print res_dir

zss = [ read_zeta_from_dir(d) for d in os.listdir(res_dir)]
zss = oc.uniq_opt_cbf(zss, 0.001)

rep_list = []
for zs in zss:
    tmp = [ [("x{0}r".format(i), z.real), ("x{0}i".format(i), -z.imag)]
            for (i, z) in zip(range(len(zs)), zs)]
    acc = []
    for s in tmp:
        acc.extend(s)
    rep_list.append(dict(acc))

def go_calc():
    cmd = 'echo "write_psi: true psi.out" >> opt_cbf.in &&'
    cmd+= '/Users/rei/src/git/opt_cbf/opt_cbf opt_cbf.in opt_cbf.out'
    return cmd

vp.run(
    root = os.path.dirname(__file__),
    files = [ "opt_cbf.in" ],
    params = vp.params.direct(rep_list),
    commands = vp.command(go_calc()),
    out_dir = vp.out_dir.flat_from("each_calc_results"),
    save_if = vp.save_if.default
)


    



