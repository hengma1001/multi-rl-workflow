import os
import sys
import shutil
import numpy as np

sys.path.append("../")
from multirl.sim.amber import AMBER_param
from multirl.sim.sim import Simulate
from multirl.sim.utils import dict_from_yaml, cal_rmsf

def sim_eval(yml_file, pdb=None, amber_bin=''):
    args = dict_from_yaml(yml_file)
    if not pdb:
        pdb = args['pdb_file']
    pdb, top = param(pdb, amber_bin=amber_bin)

    args['pdb_file'] = pdb
    args['top_file'] = top
    sim_path = sim(args)
    dcd = f"{sim_path}/output.dcd"

    rmsf = cal_rmsf(top, dcd)
    np.save(f'{sim_path}/rmsf.npy', rmsf)
    return rmsf


def param(pdb, **kwargs):
    host_dir = os.getcwd()

    # label for ligand identity
    pdb_code = os.path.basename(pdb)[:-4]

    # make the work dir
    work_dir = os.path.abspath(os.path.join(host_dir, 'input_' + pdb_code))
    os.makedirs(work_dir, exist_ok=True)

    # make a copy of pdb in the new dir
    pdb_copy = os.path.join(work_dir, os.path.basename(pdb))
    shutil.copy2(pdb, pdb_copy)

    # run and get the parameters
    os.chdir(work_dir)
    amberP = AMBER_param(pdb_copy, forcefield='ff14SB',
            watermodel='tip3p', **kwargs)
    print(amberP.prot_files, amberP.lig_files)
    pdb, top = amberP.param_comp()
    os.chdir(host_dir)
    return pdb, top

def sim(args): 
    sim = Simulate(**args)
    return sim.md_run()

