import os
from multirl.sim.evb import evb_run

if __name__ == '__main__':
    pdb = os.path.abspath('2pwz_G.pdb')
    lig_yml = os.path.abspath('lig.yml')
    ref_pdb = os.path.abspath('./5mdh_b.pdb')
    template_yml = os.path.abspath("template.yml")
    pymol_exec = '/homes/heng.ma/miniconda3/envs/pymol/bin/pymol'

    org_path = os.path.abspath('./')
    temp_dir = "evb_test"
    os.makedirs(temp_dir, exist_ok=True)
    os.chdir(temp_dir)
    md_ymls = evb_run(pdb, ref_pdb, lig_yml, template_yml, pymol_exec)
    print(md_ymls)
    os.chdir(org_path)
