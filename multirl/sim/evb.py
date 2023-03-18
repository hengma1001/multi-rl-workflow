import os
import sys
import numpy as np
import parmed as pmd
import MDAnalysis as mda

from MDAnalysis.analysis import align

from .utils import build_logger
from .utils import run_and_save
from .utils import BaseSettings
from .run import param

logger = build_logger()

class top_setup(BaseSettings):
    pdb : str
    top : str

class lig_setup(BaseSettings):
    mr_lig : top_setup
    mp_lig : top_setup


def evb_setup(pdb, ref_pdb, lig_yml, pymol_exec): 
    prot_label = os.path.basename(pdb)[:-4]
    pdb_aligned = pymol_align(pdb, ref_pdb, pymol_exec=pymol_exec)
    
    # replace his195 to hip
    prot_pdb = _mdh_mod_hip(pdb_aligned)
    # build top
    prot_crd, prot_top = param(prot_pdb, add_sol=False)
    prot_top = pmd.load_file(prot_top, xyz=prot_crd)

    # loop ligs
    ligs = lig_setup.from_yaml(lig_yml)
    for i, lig in enumerate(ligs._iter()):
        lig_type = lig[0]
        lig_top = lig[1].top
        lig_pdb = lig[1].pdb
        lig_pmd = pmd.load_file(lig_top, xyz=lig_pdb)
        comp_top = prot_top + lig_pmd
        
        if i == 0:
            sol_top = build_sol(prot_pdb, lig_pdb, lig_param_path=os.path.dirname(lig_top))
        comp_top = comb_top(comp_top, sol_top)

        output_dir = f"input_{prot_label}_{lig_type}"
        os.makedirs(output_dir)
        comp_top.save(f"{output_dir}/{prot_label}_{lig_type}.top")
        comp_top.save(f"{output_dir}/{prot_label}_{lig_type}.pdb")
        comp_top.save(f"{output_dir}/{prot_label}_{lig_type}.gro")


def _mdh_mod_hip(pdb, resid_hip = 195):
    mda_u = mda.Universe(pdb)
    res_his = mda_u.select_atoms('resname HIS')
    if len(res_his) != 0:
        hip_ind = np.argmin([np.abs(resid_hip - i.resindex) for i in res_his.residues])
        if np.abs(hip_ind - resid_hip) < 50:
            res_his.residues[hip_ind].resname = 'HIP'
            mda_u.atoms.write(pdb)
    return pdb

def pymol_align(moble, target, pymol_exec=None): 
    output = moble[:-4] + '_a.pdb'
    if not pymol_exec: 
        pymol_exec = 'pymol'
    cmd = f"{pymol_exec} -cq {os.path.dirname(__file__)}/pymol_align.py -- {moble} {target} {output}"
    with open("pymol_log", 'w') as log: 
        run_and_save(cmd, log)
    return output

def build_sol(prot_pdb, lig_pdb, lig_param_path=None): 
    prot_top = pmd.load_file(prot_pdb)
    lig_pmd = pmd.load_file(lig_pdb)
    comp_top = prot_top + lig_pmd
    
    temp_pdb = 'sol.pdb'
    comp_top.save(temp_pdb, overwrite=True)
    pdb, top = param(temp_pdb, lig_param_path=lig_param_path)
    top_sol = pmd.load_file(top, xyz=pdb)
    return top_sol
    
def comb_top(comp_top, sol_top):
    mda_comp = mda.Universe(comp_top)
    sol_comp = mda.Universe(sol_top)

    align.alignto(mda_comp, sol_comp, select='protein and name CA')
    comp_top.posistions = mda_comp.atoms.positions

    top_sol_wat = sol_top[len(mda_comp.atoms)+1:]
    comp_top += top_sol_wat
    return comp_top

if __name__ == '__main__': 
    pass