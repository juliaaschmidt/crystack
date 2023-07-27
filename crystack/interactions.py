"""
Interactions
====

    1. π-π stacking
    2. hydrophobic interactions (Car with other Car)
    3. short contacts (<= sum vdW radii)

Currently not supported:
    4. halogen interactions? to C-X...O-R or ..N-R (currently suspended)
    5. hydrogen bonds

"""
import itertools
import re
from openbabel import openbabel, pybel
from collections import namedtuple

from .tables import atomic_num, atom_vdw_radii
from .utilities import (
        euclidean3d
)
from . import config
from .stacking import PiStacks


class Dimer:

    def __init__(self, basen):
        self.basen = basen
        self.interacting_mols = re.findall(r'\d+', self.basen)
        self.omol = None
        self.all_atoms = None
        self.MolA, self.MolB = None, None
        # self.centreA, self.centreB = get_mol_centroids(self, self.interacting_mols)

    def make_mol_obj(self):
        """ creates initial Dimer object needed for analysis"""
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol", "mol")
        self.omol = openbabel.OBMol()
        obConversion.ReadFile(self.omol, f"{self.basen}.mol")
        self.all_atoms = pybel.Molecule(self.omol)
        n_atoms = self.omol.NumAtoms()
        #print("n_atoms(dimer): ", n_atoms)
        separate_mols = [m for m in self.omol.Separate()]
        self.MolA, self.MolB = pybel.Molecule(separate_mols[0]),\
                               pybel.Molecule(separate_mols[1])

        return self.omol, self.all_atoms, self.MolA, self.MolB


class Interactions(Dimer):
    """ 
    Representsl all intermolecular
    interactions in a molecular
    crystal dimer
    """
    SMALLEST_RING_SIZE = 4
    LARGEST_RING_SIZE = 6

    def __init__(self, basen):
        super().__init__(basen)
        self.short_stacks = []
        self.hydro_stacks = []
        self.pi_stacks = []

    @staticmethod
    def _get_sum_vdw_radii(atomA, atomB):
        """
        Computes the sum of the vdw radii of two atoms a and b
        a,b must be openbabel atom types
        """
        labelA = re.sub(r"\d*$", "", atomA.OBAtom.GetType().
                        replace('ar', '').replace("\d+", ""))
        labelB = re.sub(r"\d*$", "", atomB.OBAtom.GetType()
                        .replace('ar', '').replace("\d+", ""))
        vdw_dist = atom_vdw_radii[labelA] + atom_vdw_radii[labelB]
        return labelA, labelB, vdw_dist

    def short_contacts(self, vdwextra=0.0):
        """
        Computes short contacts within the dimer,
        i.e. atom-atom distance, which are lower than
        the sum of their vdW radii.

        optional arguments:
        vdwextra : if you want to include short contacts,
        with a sum larger than the sum of
        their vdW radii. e.g. vdwextra can be set to +1.0 A
        """
        self.short_stacks = []
        data = namedtuple('vdw', 'atomA idxA typeA atomB idxB typeB distance')
        atomsA = [a for a in self.MolA]
        atomsB = [a for a in self.MolB]

        for a, b in itertools.product(atomsA, atomsB):

            labelA, labelB, vdw_dist = self._get_sum_vdw_radii(a, b)
            distance = euclidean3d(a.coords, b.coords)

            if distance <= (vdw_dist + vdwextra):
                distance_diff = round(distance-(vdw_dist + vdwextra), 3)

                short_contact = data(atomA=a, idxA=a.idx, typeA=labelA,
                                     atomB=b, idxB=b.idx, typeB=labelB,
                                     distance=distance)

                self.short_stacks.append(short_contact)

                print(f"{a.idx} {labelA} <-> {b.idx} {labelB} \
                d:{distance} Lower by:{distance_diff} \
                vdWthresholdIncrease:{vdwextra}")
        return self.short_stacks

    @staticmethod
    def _get_hydrophobic_atoms(molecule):
        """
        Select all carbon atoms which have only carbons
        and/or hydrogens connected to them.

        Returns:
        hydro_candidates : set of candidate atoms, capable
        of hydrophobic interactions.
        """
        data = namedtuple('hydrophobic', 'atom type idx')
        hydro_atm = [a for a in molecule if a.atomicnum == int(atomic_num['C'])
                     and set([natom.GetAtomicNum() for natom
                     in pybel.ob.OBAtomAtomIter(a.OBAtom)])
                     .issubset({atomic_num['H'], atomic_num['C']})]

        hydro_candidates = [data(atom=atom, type=atom.OBAtom.GetType(),
                            idx=atom.idx) for atom in hydro_atm]
        return hydro_candidates

    def hydrophobic_interactions(self):
        """
        Calculate hydrophobic interactions between molecules

        Returns
        -------
        :class:`list`
            list of namedtuples with all
            stacking information

        References
        ----------
        Definition: All pairs of qualified carbon atoms within a
        distance of HYDROPH_DIST_MAX taken from plip module:
        https://github.com/pharmai/plip
        """
        data = namedtuple('hydroph_interaction',
               'atomA idxA typeA atomB idxB typeB distance')
        hydro_setA = self._get_hydrophobic_atoms(self.MolA)
        hydro_setB = self._get_hydrophobic_atoms(self.MolB)

        for a, b in itertools.product(hydro_setA, hydro_setB):

            e = euclidean3d(a.atom.coords, b.atom.coords)

            if config.MIN_DIST < e < config.HYDROPH_DIST_MAX:
                contact = data(atomA=a.atom, idxA=a.idx, typeA=a.type,
                               atomB=b.atom, idxB=b.idx, typeB=b.type,
                               distance=e)
                self.hydro_stacks.append(contact)

                print(f"Hydro. int: {a.idx} {a.type}-{b.idx} {b.type} {e}")
        return self.hydro_stacks

    def calculate_hydrophobic_interactions(self):
        """
        Calculate hydrophobic interactions between
        neighbouring molecules

        Returns
        -------
        :class:`list`
            list of namedtuples with all
            stacking information
        """
        return self.hydrophobic_interactions()

    def calculate_short_contacts(self, **kwargs):
        """
        Calculate short contacts between molecules

        Returns
        -------
        :class:`list`
            list of namedtuples with all
            short contacts
        """
        return self.short_contacts(**kwargs)

    def calculate_pi_stacks(self):
        """
        Calculate π-stacks between molecules

        Returns
        -------
        :class:`list`
            list of namedtuples with all
            stacking information
        """
        stacks = PiStacks(self.omol, self.all_atoms)
        stacks.find_rings()
        return stacks.pistacking()

    def calculate_all_contacts(self, vdwextra=0.):
        """ 
        Executes all interaction functionality:

        1. π-π stacks
        2. short contacts
        3. hydrophobic interactions
        (4.) halogen interactions (currently suspended)
        """
        # calculate π-stacks
        stacks = PiStacks(self.omol, self.all_atoms)
        stacks.find_rings()
        self.pi_stacks = stacks.pistacking()
        # stacks.create_dummy_centroids(self.interacting_mols)

        self.hydrophobic_interactions()
        self.short_contacts(vdwextra=0.0)
        return self.pi_stacks, self.hydro_stacks, self.short_stacks

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'{self.__class__.__name__}(pi_stacks={self.pi_stacks})'