import glob
import numpy as np
from multirl.sim.evb import evb_analysis

if __name__ == '__main__':
    md_paths = glob.glob("/lambda_stor/homes/heng.ma/Research/gene_transformer/mdh_evb/md_runs/run_ymls/evb_32a/md_run/md_run_*")
    wham = '/lambda_stor/homes/heng.ma/Research/gene_transformer/mdh/evb/EVB-amber/examples/wham/wham/wham'
    pmf = evb_analysis(
        md_paths,
        wham_exe=wham)
    np.save('pmf.npy', pmf)
