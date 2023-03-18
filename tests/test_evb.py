import sys
sys.path.append("../")
from multirl.sim.evb import evb_setup

if __name__ == '__main__':
    pdb = '2pwz_G.pdb'
    lig_yml = 'lig.yml'
    ref_pdb = './5mdh_b.pdb'
    pymol_exec = '/homes/heng.ma/miniconda3/envs/pymol/bin/pymol'
    evb = evb_setup(pdb, ref_pdb, lig_yml, pymol_exec)