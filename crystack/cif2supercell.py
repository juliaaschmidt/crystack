from re import I
import numpy as np
import sys
import pandas as pd
from scipy.spatial import distance
import subprocess as sp

import pywindow as pw
from pymatgen.io.cif import CifParser
from pymatgen.io.xyz import XYZ
from pymatgen.core.structure import Structure

from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import JmolNN  # used for deciding which atoms are bonded
#from utilities import centroid
from .utilities import euclidean3d
import csv
"""

PYMATGEN ALTERNATIVE TO BUILDING THE SUPERCELL FROM CIF FILE

requirements:
* ! only works in the py36 env
- pymatgen
- pywindow (python 3.6, git clone from repos)

ideas to improve this:
- get the neighbouring dimers according to axis directions (?)

INPUT:
- cif file
- n_atoms (single molecule)

OUTPUT:
- dimer xyz (nonreal)
"""
def make_header(n_atoms):
    dimer_n_atoms = str(n_atoms * 2)
    with open('header_dimer', 'w') as the_file:
        the_file.write(f"{dimer_n_atoms}\n\n")

# def convert_xyz_mol():

class Crystal:

    def __init__(self, filn, n_atoms):
        self.filn = filn
        self.basen = filn.split(".")[0]
        self.n_atoms = n_atoms

    def create_supercell(self, size=[4,4,4]):
        """
        Creates a n x n x n supercell based on
        the cif file.

        Parameters
        -------
        supercell_size : :class:`list` (optional)
            3-D vector specifying the size of the supercell,
            default is [4,4,4] i.e. 4x4x4 supercell

        Returns
        -------
        :class: `xyz file`
            file with cartesian coordinates of all the
            atoms in the supercell
        """
        scaling_matrix = [
            [size[0], 0, 0],
            [0, size[1], 0],
            [0, 0, size[2]]
        ]
        s = Structure.from_file(self.filn)
        s.make_supercell(scaling_matrix, to_unit_cell=False)
        self.xyzrep = XYZ(s) # convert to xyz
        self.basen = self.filn.split(".")[0]
        self.xyzrep.write_file(f"{self.basen}_supercell.xyz")  # write supercell to file
        # print(np.shape(self.xyzrep.molecule))
        # print(self.xyzrep.molecule)

    @staticmethod
    def distance(self):
        return np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)  

    def find_molecules(self):
        #         Rcov(x) + Rcov(y) - t < R(x,y) < Rcov(x) + Rcov(y) + t
        # where Rcov is the covalent radius and the tolarenace (t) is set to 0.4 Angstrom.
        n, label, x, y, z = self.xyzrep.split
        array = np.array(x,y,z)
        dist=[[0]*len(array[0])]*len(array)  
        for i in range(len(array)-1):  
            for j in range(i+1,len(array)):  
                dist[i][j]=distance(array[i],array[j])  

        # d = ((x2 - x1)2 + (y2 - y1)2 + (z2 - z1)2)1/2   
        # n, label, x, y, z = self.xyzrep
        pass

    # # find centre-of-mass of supercell and return all COMs and mols xyz.
    # centroid, coms_list, all_mols = find_supercell_centre(basen)

    def find_supercell_centre(self):
        """
        Takes the supercell (XYZ) and finds the centre-of-mass (COM)
        of the supercell.

        Returns
        -------
        :class:`list`
            3D coordinates of the global centroid of the supercell
        :class:`list`
            list of all centre-of-mass (COMs) computed for each not-broken
            molecule in the supercell
        :class: `list`
            all_mols: list of xyz of all molecules in supercell.

        !TODO: This is the only method requiring pywindow. Should be equally
        !feasible with pymatgen.
        """
        molsys = pw.MolecularSystem.load_file(f"{self.basen}_supercell.xyz")
        molsys.make_modular()  # finds separate molecules in supercell

        center_center = []
        # list of all the COMs of the mols in supercell (locations.txt)
        all_mols = []
        for molecule in molsys.molecules:
            mol = molsys.molecules[molecule]
            if mol.no_of_atoms == self.n_atoms:  # check if mol is broken
                info = mol.calculate_centre_of_mass()
                center_center.append(info.tolist())  # com to locations list
                all_mols.append(mol)  # molecule to mols list

        center_center = list(center_center)
        n_mols = len(center_center)
        print(f"There are {n_mols} molecules with the correct #atoms in supercell.")
        assert n_mols != 0, "Error: Incorrect n_atoms per mol entered"
        centroid = np.mean(np.array(center_center), axis=0)
        return centroid, center_center, all_mols

    @staticmethod
    def closest_node(node, nodes):
        """
        Finds the molecule closest to the centre of mass of
        the supercell (closest distance by COM-COM distance).

        Returns
        -------
        :class:`list`
            position of molecule closest to com of supercell
        :class:`int`
            idx of molecule closest to com of supercell
        """
        closest_index = distance.cdist([node], nodes).argmin()
        return nodes[closest_index], closest_index


def duplicate_dimers_finder(dist_only):
    """
    Checks whether current dimer distance
    is same as previous distance.

    Returns
    -------
    :class:`list`
        list of indices of all duplicate dimer distances
    * TODO: there are easier way to do this in pandas, like keepfirst,
            * but didn't work.
    """
    dist_only = dist_only.round(6)  # avoid differences due to numerical noise
    idx_duplicates = []
    last_item = None

    for i, dist in dist_only.iteritems():
        if dist == last_item:
            idx_duplicates.append(i)
            last_item = dist  # needs to be inside if statement!
    return idx_duplicates


def find_closest_packed_dimers(closest, closest_idx, coms_list):
    """

    We only require the closest dimers, not dimers which are further apart
    so in reality they would have another molecule sitting in between.

    Hence we sort all mols by their distance from our most central molecule
    Then, we only select the 5 closest unique dimers, ignoring dimers that
    are duplicates ie. have the same COM-COM distance.

    returns
        - 5 unique nearest neighbour dimers
    """
    # find all mols closest to the one closest to the centroid of the supercell
    df = pd.DataFrame(np.array(coms_list), columns=['x', 'y', 'z'])
    df["dist"] = np.sqrt(((df.x-closest[0])**2) +
                        ((df.y-closest[1])**2)+((df.z-closest[2])**2))

    # idx of df indicates molecules to be used/at closest distance
    sorted_df = df.sort_values(by=["dist"])

    dist_only = pd.Series(sorted_df['dist'], dtype="float64")
    idx_duplicates = duplicate_dimers_finder(dist_only)

    # print(sorted_df)
    # remove duplicate dimer motifs with the same COM-COM distance:
    clean = sorted_df.drop(idx_duplicates)
    # clean = sorted
    # print(idx_duplicates)
    # save unique positions to file (optional)
    sorted_df.to_csv("coms_distance_sorted.csv")
    return clean


def generate_sm_xyz(central_top, all_mols):
    # Generates single molecule real_xyz files
    for idx in central_top:
        print(all_mols[idx])
        all_mols[idx].dump_molecule(
            "./mol_{0}.pdb".format(idx),
            include_coms=False,
            override=True)
        sp.call(f"obabel -ipdb mol_{idx}.pdb -oxyz mol_{idx}.xyz > mol_{idx}.xyz", shell=True)
    sp.call("rm mol_*.pdb", shell=True)


def write_dimer_xyz():
    with open("file.xyz", 'w') as xyz_file:
        xyz_file.write("%d\n%s\n" % (len(molecule.atoms), title))
        for atom in molecule.atoms:
            xyz_file.write("{:4} {:11.6f} {:11.6f} {:11.6f}\n".format(
                atom.atomic_symbol, atom.coordinates.x, atom.coordinates.y, atom.coordinates.z))

def generate_dimer_xyz(idx_closest, central_top):
    """
    Function generates real dimer xyz files and opens them.
    """
    # remove header from central molecule file
    sp.call(
        f"tail -n +3 mol_{idx_closest}.xyz > nonreal_{idx_closest}.xyz",
        shell=True)

    # makes dimer files from monomers
    orient_dict = {}
    for idx, name in enumerate(central_top[1:]):  # start from 2nd mol in list
        sp.call(f"tail -n +3 mol_{name}.xyz > nonreal_{name}.xyz", shell=True)
        # rm header

        # make nonreal dimer xyz
        sp.call(
            f"cat nonreal_{idx_closest}.xyz nonreal_{name}.xyz > \
                dimer_{idx_closest}_{name}.xyz", shell=True)

        #  make real dimer xyz
        sp.call(
            f"cat header_dimer dimer_{idx_closest}_{name}.xyz > \
                real_dimer_{idx_closest}_{name}.xyz", shell=True)
        
        dominant_dir = find_abc_orientation(name, idx_closest)
        orient_dict[name] = dominant_dir

        sp.call(f"obabel -ixyz real_dimer_{idx_closest}_{name}.xyz -omol real_dimer_{idx_closest}_{name}.mol > real_dimer_{idx_closest}_{name}.mol", shell=True)

        sp.call(f"rm real_dimer_{idx_closest}_{name}.xyz", shell=True)
        sp.call(f"rm dimer_{idx_closest}_{name}.xyz", shell=True)
        sp.call(f"rm nonreal_{name}.xyz", shell=True)
        sp.call(f"rm mol_{name}.xyz", shell=True)

    return orient_dict

def run_dimer_generator(filn, n_atoms, out_mol=5):
    # OUTPUT settings:
    # out_mol = 5
    # number of unique nearest neighbour dimers you want to
    # extract from crystal

    # Input handling
    #filn = sys.argv[1]
    #self.basen = filn.split(".")[0]
    #ibasen = filn.split(".")[0]
    #n_atoms = int(sys.argv[2]) # n atoms per mol

    make_header(n_atoms)
    crystal_to_generate = Crystal(filn, n_atoms)
    crystal_to_generate.create_supercell()

    # crystal_to_generate.find_molecules()

    # find centre-of-mass of supercell and return all COMs and mols xyz.
    supercell_centroid, coms_list, all_mols = crystal_to_generate.find_supercell_centre()

    # find molecule closest to centroid of supercell
    closest, idx_closest = crystal_to_generate.closest_node(supercell_centroid, coms_list)

    # get list of all the unique closest dimers in ascending order
    clean = find_closest_packed_dimers(closest, idx_closest, coms_list)

    # make dimer xyz files for (25 or) 15 closest molecules
    central_top = clean[:out_mol+1].index.values.tolist()

    # makes single molecule xyz files (rm pdb files created during this)
    generate_sm_xyz(central_top, all_mols)

    # make dimer xyz
    orient_dict = generate_dimer_xyz(idx_closest, central_top)
    # # use the closest 5-6 mols as dimers
    # with open('dominant_orientations.csv', 'w') as csv_file:  
    # writer = csv.writer(csv_file)
    # for key, value in orientations_dict.items():
    #     writer.writerow([key, value])

    sp.call("rm header_dimer", shell=True)
    return orient_dict


def find_abc_orientation(idx1, idx_origin):
    coms_file = pd.read_csv("coms_distance_sorted.csv", index_col=0)
    # print(coms_file.head())
    row1 = coms_file.loc[idx1]
    row2 = coms_file.loc[idx_origin]
    # print(row1)
    # print(row2)

    orientation_vec = [row1.x- row2.x, row1.y - row2.y, row1.z - row2.z]
    # print(idx1, idx_origin)

    absolute_orientations = list(map(abs, orientation_vec))
    max_value = max(absolute_orientations)
    max_index = absolute_orientations.index(max_value)
    #print(max_index)

    if max_index == 0:
        # print("a")
        return "a"
    elif max_index == 1:
        # print("b")
        return "b"
    else: 
        # print("c")
        return "c"
    
