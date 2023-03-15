import sys
sys.path.append("../")
from multirl.sim.run import main as sim_eval

if __name__ == '__main__':
    pdb = '2pwz_G.pdb'
    yml = 'md.yml'
    rmsf = sim_eval(yml, pdb=pdb)
    print(rmsf)