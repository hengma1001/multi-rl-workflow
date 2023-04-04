import os
import tempfile
import subprocess
import numpy as np
import pandas as pd
import parmed as pmd
import MDAnalysis as mda

from MDAnalysis.analysis import align

from .utils import build_logger
from .utils import BaseSettings
from .utils import run_and_save, create_path
from .utils import dict_from_yaml, dict_to_yaml
from .utils import run_in_tempdir
from .run import param, sim

logger = build_logger()


class top_setup(BaseSettings):
    pdb: str
    top: str


class lig_setup(BaseSettings):
    mr_lig: top_setup
    mp_lig: top_setup


def evb_run(pdb, ref_pdb, lig_yml, template_yml,
            pymol_exec='pymol', wham_exe='wham_exe',
            amber_bin=''):
    """Run empirical valence bond method for reaction free energy profile

    Args: 
        pdb: Initial structure
        ref_pdb: pdb pose with binding ligands 
        lig_yml: ligand top setup for both reactant and product
        template_yml: template setup for evb sim runs 
        pymol_exec: pymol executable path
        wham_exec: wham executable path
        amber_bin: amber bin directory
    Returns:

    """
    logger.info("Building pdb and top...")
    sim_setups = evb_setup(pdb, ref_pdb, lig_yml, pymol_exec, amber_bin)
    # build yml files for md simulations
    logger.info("Building MD setup for each sampling point...")
    md_ymls = evb_ymls(template_yml, sim_setups)

    # return md_ymls
    # TODO: added function to run and analyze simulations
    sim_paths = []
    for yml in md_ymls:
        sim_paths.append(sim(**dict_from_yaml(yml)))

    # gather and analysis
    pmf = evb_analysis(
        sim_paths, wham_exe=wham_exe
    )
    return pmf


def evb_ymls(template_yml, sim_setups):
    run_setup = dict_from_yaml(template_yml)
    for sim_name in sim_setups:
        sim_sys = sim_name.split('_')[0]
        label = os.path.basename(sim_setups[sim_name])[:-3]
        run_setup['output_dir'] = f"evb_{label}"
        for filetype in ["top", "pdb"]:
            run_setup['md_setup'][f'{sim_sys}_{filetype}'] = \
                os.path.abspath(f"{sim_setups[sim_name]}.{filetype}")

    dict_to_yaml(run_setup, 'test.yml')
    return build_ymls(run_setup)


def evb_setup(pdb, ref_pdb, lig_yml, pymol_exec, amber_bin):
    prot_label = os.path.basename(pdb)[:-4]
    pdb_aligned = pymol_align(pdb, ref_pdb, pymol_exec=pymol_exec)

    # replace his195 to hip
    prot_pdb = _mdh_mod_hip(pdb_aligned)
    # build top
    prot_crd, prot_top = param(prot_pdb, add_sol=False, amber_bin=amber_bin)
    prot_top = pmd.load_file(prot_top, xyz=prot_crd)

    # loop ligs
    ligs = lig_setup.from_yaml(lig_yml)
    sim_setups = {}
    for i, lig in enumerate(ligs._iter()):
        lig_type = lig[0].replace('_lig', '')
        lig_top = lig[1].top
        lig_pdb = lig[1].pdb
        lig_pmd = pmd.load_file(lig_top, xyz=lig_pdb)
        comp_top = prot_top + lig_pmd

        if i == 0:
            sol_top = build_sol(prot_pdb, lig_pdb, lig_param_path=os.path.dirname(lig_top), 
                    amber_bin=amber_bin)
        comp_top = comb_top(comp_top, sol_top)

        output_dir = f"input_{prot_label}_{lig_type}"
        os.makedirs(output_dir)

        filesave = f"{output_dir}/{prot_label}_{lig_type}"
        comp_top.save(f"{filesave}.top")
        comp_top.save(f"{filesave}.pdb")
        comp_top.save(f"{filesave}.gro")
        sim_setups[lig_type] = filesave
    logger.info(f"Finished setup for {filesave}")
    return sim_setups


def _mdh_mod_hip(pdb, resid_hip=195):
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


def build_sol(prot_pdb, lig_pdb, lig_param_path=None, **kwargs):
    """Add water to the prot+lig complex"""
    prot_top = pmd.load_file(prot_pdb)
    lig_pmd = pmd.load_file(lig_pdb)
    comp_top = prot_top + lig_pmd

    host_dir = os.getcwd()
    try:
        with tempfile.TemporaryDirectory() as tp:
            os.chdir(tp)
            temp_pdb = f'{tp}/sol.pdb'
            comp_top.save(temp_pdb, overwrite=True)
            pdb, top = param(temp_pdb, lig_param_path=lig_param_path, **kwargs)
            top_sol = pmd.load_file(top, xyz=pdb)
    finally:
        os.chdir(host_dir)
    return top_sol


def comb_top(comp_top, sol_top):
    mda_comp = mda.Universe(comp_top)
    sol_comp = mda.Universe(sol_top)
    align.alignto(mda_comp, sol_comp, select='protein and name CA')
    comp_top.positions = mda_comp.atoms.positions

    top_sol_wat = sol_top[len(mda_comp.atoms)+1:]
    comp_top += top_sol_wat
    comp_top.box = sol_top.get_box()

    return comp_top


def build_ymls(evb_setup):
    evb_cfg = evb_setup['evb_cfg']
    rc_min, rc_max = evb_cfg['rc_min'], evb_cfg['rc_max']
    rc_inc = evb_cfg['rc_inc']
    rc0_list = np.arange(rc_min, rc_max+rc_inc, rc_inc)

    md_setup = evb_setup['md_setup']
    # correcting file path
    input_files = ['mr_pdb', 'mp_pdb', 'mp_top', 'mr_top']
    for input in input_files:
        if input in md_setup and md_setup[input]:
            if not os.path.isabs(md_setup[input]):
                md_setup[input] = \
                    os.path.join(os.getcwd(), md_setup[input])
                logger.debug(f"updated entry{input} to {md_setup[input]}.")
    md_path = create_path(dir_type='md')

    md_ymls = []
    for rc0 in rc0_list:
        for run in ['mr', 'mp']:
            md_yml = f"{md_path}/md_{rc0:.5f}_{run}.yml"
            md_setup_copy = md_setup.copy()
            md_setup_copy['pdb_file'] = md_setup_copy[f'{run}_pdb']
            md_setup_copy['top_file'] = md_setup_copy[f'{run}_top']
            # protein atom number
            mda_u = mda.Universe(md_setup_copy['pdb_file'])
            protein = mda_u.select_atoms('protein')
            n_atoms_protein = protein.n_atoms
            # umb setup
            dbonds_umb = evb_cfg['dbonds_umb'].copy()
            dbonds_umb['rc0'] = float(rc0)
            dbonds_umb['atom_i'] = int(evb_cfg[f"{run}_setup"]['mr_atom']) - 1 + n_atoms_protein
            dbonds_umb['atom_j'] = int(evb_cfg[f"{run}_setup"]['mp_atom']) - 1 + n_atoms_protein
            dbonds_umb['atom_k'] = int(evb_cfg[f"{run}_setup"]['h_atom']) - 1 + n_atoms_protein
            md_setup_copy['dbonds_umb'] = dbonds_umb
            # morse setup
            morse_bond = evb_cfg['morse_bond'].copy()
            morse_bond['atom_i'] = int(evb_cfg[f"{run}_setup"][f'{run}_atom']) - 1 + n_atoms_protein
            morse_bond['atom_j'] = int(evb_cfg[f"{run}_setup"]['h_atom']) - 1 + n_atoms_protein
            md_setup_copy['morse_bond'] = morse_bond
            dict_to_yaml(md_setup_copy, md_yml)
            md_ymls.append(md_yml)

    return md_ymls


def evb_analysis(
        md_paths: list,
        skip_start: int = 100,
        wham_exe: str = 'wham',):
    evb_df = []
    for md in md_paths:
        rc_0 = float(os.path.basename(md).split('_')[-1])
        md_type = os.path.basename(md).split('_')[-2]

        # sim_log = f"{md}/output.log"
        # sim_df = pd.read_csv(sim_log)
        rc_log = f"{md}/output.rc"
        rc_df = pd.read_csv(rc_log)

        local_df = {
            "rc0": rc_0,
            "rc": rc_df['rc'].to_list()[skip_start:],
            "dist_mr": rc_df['dist_mr'].to_list()[skip_start:],
            "dist_mp": rc_df[' dist_mp'].to_list()[skip_start:],
            'frame': list(np.arange(len(rc_df)))[skip_start:],
            # 'E_p': sim_df['Potential Energy (kJ/mole)'].to_list()[skip_start:],
            'run_type': md_type}
        evb_df.append(local_df)

    # dataframe
    df = pd.DataFrame(evb_df)
    df = df.explode(column=["rc", "dist_mr", "dist_mp", "frame"]).reset_index(drop=True)
    # df = df.explode(column=['rc', "dist_mr", "dist_mp", 'frame', 'E_p']).reset_index(drop=True)
    df = df.astype({'rc0': float, 'rc': float, 'frame': int,})
    df['E_umb'] = 1/2 * 500000 * (df.rc/10 - df.rc0)**2

    rmsf = run_wham(df, wham_exe=wham_exe)
    return rmsf


@run_in_tempdir
def run_wham(df, wham_exe):
    data_path = 'data'
    os.makedirs(data_path, exist_ok=1)
    with open('meta.dat', 'w') as f_meta:
        for i, rc0 in enumerate(df.rc0.unique()):
            sub_df = df[df.rc0 == rc0][['frame', 'rc']]
            save_sim = f"{data_path}/run_{i}"
            np.savetxt(save_sim, sub_df.to_numpy())
            f_meta.write(f"{save_sim} {rc0*10:.2f} {5000}\n")

    wham_cmd = f'{wham_exe} -.5 .5 50 0.0000000001 300 0 meta.dat pmf_out 100 42'
    subprocess.check_output(wham_cmd, shell=True)
    pmf = np.loadtxt('pmf_out')

    return pmf


if __name__ == '__main__':
    pass
