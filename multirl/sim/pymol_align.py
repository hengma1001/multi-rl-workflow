import os
import sys
import pymol

# print(sys.argv)
pdb, pdb_ref, pdb_out = sys.argv[1:]

pymol.cmd.load(pdb_ref, 'ref')
pdb_label = os.path.basename(pdb)
pymol.cmd.load(pdb, pdb_label)

pymol.cmd.align(pdb_label, 'ref')
pymol.cmd.save(pdb_out, pdb_label)