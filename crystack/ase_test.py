from ase import neighborlist
from ase.build import molecule
from scipy import sparse

import sys
import pandas as pd
from scipy.spatial import distance
import subprocess as sp

import pywindow as pw
from pymatgen.io.cif import CifParser
from pymatgen.io.xyz import XYZ
from pymatgen.core.structure import Structure
from pymatgen.core.structure import extract_cluster

# Input handling
filn = sys.argv[1]
basen = filn.split(".")[0]
n_atoms = int(sys.argv[2])  # make sure to input correct thing here
# parser = CifParser(filn)
# struc = parser.get_structures()[0]

s = Structure.from_file(filn)
scaling_matrix = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]  # [[4, 0, 0], [0, 4, 0], [0, 0, 4]]  # create 4x4x4 matrixs.make_supercell(scaling_matrix, to_unit_cell=False)

s.make_supercell(scaling_matrix, to_unit_cell=False)
extract_cluster(s)
print(s)
cutOff = neighborlist.natural_cutoffs()
neighborList = neighborlist.NeighborList(cutOff, self_interaction=False, bothways=True)
neighborList.update(s)
matrix = neighborList.get_connectivity_matrix()
#or: matrix = neighborlist.get_connectivity_matrix(neighborList.nl)
n_components, component_list = sparse.csgraph.connected_components(matrix)
idx = 1
molIdx = component_list[idx]
print("There are {} molecules in the system".format(n_components))
print("Atom {} is part of molecule {}".format(idx, molIdx))
molIdxs = [ i for i in range(len(component_list)) if component_list[i] == molIdx ]
print("The following atoms are part of molecule {}: {}".format(molIdx, molIdxs))