import os 
import time
import yaml
import logging
import argparse
import tempfile
import subprocess
import warnings

from typing import Union
from pathlib import Path
from typing import Type, TypeVar
# from pydantic import BaseSettings as _BaseSettings

import parmed as pmd
import MDAnalysis as mda

from rdkit import Chem
from mendeleev import element
from MDAnalysis.analysis import align, rms


PathLike = Union[str, Path]
_T = TypeVar("_T")


def build_logger(debug=0):
    logger_level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(level=logger_level, format='%(asctime)s %(message)s')
    logger = logging.getLogger(__name__)
    return logger

def dict_from_yaml(yml_file): 
    return yaml.safe_load(open(yml_file, 'r'))

def dict_to_yaml(dict_t, yml_file): 
    with open(yml_file, 'w') as fp: 
        yaml.dump(dict_t, fp, default_flow_style=False)

class yml_base(object): 
    def dump_yaml(self, cfg_path: PathLike) -> None: 
        dict_to_yaml(self.get_setup(), cfg_path)

def create_path(dir_type='md', sys_label=None, create_path=True): 
    """
    create MD simulation path based on its label (int), 
    and automatically update label if path exists. 
    """
    dir_path = f'{dir_type}_run'
    if sys_label: 
        dir_path = f'{dir_path}_{sys_label}'
    if create_path:
        os.makedirs(dir_path, exist_ok=True)
    return os.path.abspath(dir_path)

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--config", help="YAML config file", type=str, required=True,
    )
    parser.add_argument(
        "-d", "--dry_run", help="Option to only build all the ymls", 
        type=bool, default=False, required=False,
    )
    args = parser.parse_args()
    return args


def trim_line(line):
    return line.split()[0].replace('"', '')


def find_diff(a, b):
    """find elements that are in A, but not in B
    """
    return sorted([i for i in a if i not in b])


def run_and_save(command, log):
    tsk = subprocess.Popen(
            command,
            stdout=log, # subprocess.PIPE,
            stderr=subprocess.STDOUT,
            shell=True)
    tsk.wait()


def get_mismatch(a, b):
    return find_diff(a, b), find_diff(b, a)


def atomList_almostEqual(a, b):
    if len(a) == len(b):
        a, b = get_atomname(a), get_atomname(b)
        # b = sorted([i[:-1] if len(i) > 1 else i for i in b])
        x = [i != j for i, j in zip(a, b)]
        if sum(x) == 0:
            return True
    return False


def get_atomname(a):
    """remove the last label in atomnames"""
    return sorted([i[:-1] if len(i) > 1 else i for i in a])


def get_ligand(pdb_file):
    mda_trj = mda.Universe(pdb_file)
    lig = mda_trj.select_atoms('not protein')
    pdb_lig = os.path.abspath('./lig.pdb')
    lig.write(pdb_lig)
    return pdb_lig


def get_protein(pdb_file):
    mda_trj = mda.Universe(pdb_file)
    lig = mda_trj.select_atoms('protein')
    pdb_lig = os.path.abspath('./prot.pdb')
    lig.write(pdb_lig)
    return pdb_lig


def get_lig_name(lig_pdb):
    mda_u = mda.Universe(lig_pdb)
    return mda_u.atoms[0].resname

def update_pdb_obabel(pdb_file, format='pdb'):
    """
    add correct conect info to pdb structure 
    obabel -ipdb lig.pdb -opdb >  lig_obabel.pdb
    """
    pdb_ob = pdb_file[:-4] + f'_ob.{format}'
    subprocess.check_output(
        f'obabel -ipdb {pdb_file} -o{format} >  {pdb_ob}',
        shell=True)
    return pdb_ob


def get_formal_charge(pdb_file, format='pdb'): 
    pdb_file = update_pdb_obabel(pdb_file, format=format)
    if format == 'pdb':
        mol = Chem.MolFromPDBFile(pdb_file)
    elif format == 'mol2': 
        mol = Chem.MolFromMol2File(pdb_file)
    else: 
        raise Exception("Unknown format...")
    return Chem.GetFormalCharge(mol)


def get_lig_charge(pdb_file): 
    try:
        lig_charge = get_formal_charge(pdb_file, format='mol2')
    except: 
        lig_charge = get_formal_charge(pdb_file, format='pdb')
    n_electron = get_n_electron(pdb_file)
    if (n_electron % 2 == 0) & (lig_charge % 2 == 0): 
        return lig_charge 
    elif (n_electron % 2 != 0) & (lig_charge % 2 != 0):
        return lig_charge
    elif n_electron % 2 == 0: 
        warnings.warn(f"Using 0 for ligand charge. However, number of electron "\
            f"{n_electron} and charge {lig_charge} "\
            f"are mismatch for ligand {os.path.abspath(pdb_file)}")
        return 0
    else:
        raise Exception(f"Number of electron {n_electron} and charge {lig_charge} "\
            f"are mismatch for ligand {os.path.abspath(pdb_file)}")


def get_n_electron(pdb_file): 
    mda_u = mda.Universe(pdb_file)
    n_ele = [element(atom.element).atomic_number for atom in mda_u.atoms]
    return sum(n_ele)


def is_protein(pdb_file):
    mda_trj = mda.Universe(pdb_file)
    not_prot = mda_trj.select_atoms('not protein')
    if not_prot.n_atoms == 0:
        return True
    else:
        return False


def run_at_temp(func):
    """
    Run functions at a temp dir
    """
    def wrapper(*args, **kwargs):
        current_dir = os.getcwd()
        temp_path = tempfile.TemporaryDirectory()
        os.chdir(temp_path.name)
        output = func(*args, **kwargs)
        os.chdir(current_dir)
        return output
    return wrapper


def clean_pdb(pdb_file):
    """
    Remove all entris in pdb files other than `ATOM` and HETATM`
    """
    with open(pdb_file, 'r') as pdb:
        pdb_atoms = [
            line for line in pdb
            if line.startswith('ATOM') or line.startswith('HETATM')]
    with open(pdb_file, 'w') as pdb:
        pdb.write(''.join(pdb_atoms))


def to_pdb(pos_file, top_file, pdb_file):
    top = pmd.load_file(top_file, xyz=pos_file)
    top.write_pdb(pdb_file)


def missing_hydrogen(pdb_file):
    """
    Check whether a pdb file contains H atoms

    Parameters
    ----------
    pdb_file : str
        path to input pdb file

    Returns
    -------
    missingH : bool
        True if missing H, false otherwise
    """
    trj = mda.Universe(pdb_file)
    hydrogens = trj.select_atoms('name H*')
    missingH = True if hydrogens.n_atoms == 0 else False
    return missingH


def remove_hydrogen(pdb_file, pdb_noH_file):
    """
    remove H atoms from a pdb file

    Parameters
    ----------
    pdb_file : str
        path to input pdb file
    pdb_noH_file : str
        path to write pdb file with H removed
    """
    trj = mda.Universe(pdb_file)
    trj_noH = trj.select_atoms('not name H* and not name h*')
    trj_noH.write(pdb_noH_file)


def add_hydrogen(pdb_file):
    """
    add hydrogens to pdb structure if missing hydrogen atoms
    obabel -ipdb adp.pdb -h -opdb >  adph.pdb
    """
    if not missing_hydrogen(pdb_file):
        return pdb_file
    else:
        pdb_noH = pdb_file

    pdb_h = pdb_file[:-4] + '_h.pdb'
    subprocess.check_output(
        f'obabel -ipdb {pdb_noH} -h -opdb >  {pdb_h}',
        shell=True)
    clean_pdb(pdb_h)
    return pdb_h


def align_to_template(pdb_file, ref_file, pdb_output):
    """
    align frame to target
    """
    pdb = mda.Universe(pdb_file)
    ref = mda.Universe(ref_file)
    _ = align.alignto(pdb, ref, select='protein and name CA')
    pdb.atoms.write(pdb_output)


def cal_rmsf(pdb, dcd): 
    mda_u = mda.Universe(pdb, dcd)
    selection = 'protein and name CA'
    
    average_frame = align.AverageStructure(mda_u, mda_u, 
                            select=selection, ref_frame=0).run()
    ref = average_frame.results.universe

    aligner = align.AlignTraj(mda_u, ref, select=selection, in_memory=True).run()

    calphas = mda_u.select_atoms(selection)
    rmsfer = rms.RMSF(calphas, verbose=True).run()

    return rmsfer.rmsf
    
