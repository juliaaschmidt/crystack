"""
PiStacks
====
"""
import numpy as np
import itertools
from openbabel import pybel
from collections import namedtuple

from .io_tools import save_files
from .utilities import (
        euclidean3d, centroid, projection,
        vecangle, unit_vector, vector
)
from . import config


class PiStacks:
    """
    Represents all π-π stacking in a dimer object.
    """
    SMALLEST_RING_SIZE = 4
    LARGEST_RING_SIZE = 6

    def __init__(self, omol, all_atoms ):
        self.omol = omol
        self.all_atoms = all_atoms
        self.all_ring_atoms = []
        self.rings = []
        self.ring_centres = []
        self.all_ring_atoms = []
        self.pi_stacks = []

        """
        Initialise an :class:`.PiStacks` instance.

        Parameters
        ----------
        !TODO missing parameters here
        rings : :class:`list`
            A list of namedtuples with all ring
            information extracted i.e.
            atomic coordinates, normal vector
            of ring plane, centroid of ring
            plane, ring size
        ring_centres : :class:`list`
            The position of all ring centroids.
        all_ring_atoms :
        pi_stacks :
        """
    @property
    def rings_per_mol(self):
        return int(len(self.ring_centres)/2)

    @staticmethod
    def _is_ring_planar(ring, r_atoms):
        """ 
        Determine whether the ring is
        sufficiently planar to be considered
        aromatic.

        Parameters
        ----------


        References
        ----------
        Definition: aromatic ring definition
        adapted from receptor-ligand interaction
        package as presented in PLIP:
        https://github.com/pharmai/plip
        """
        normals = []
        for a in r_atoms:
            adj = pybel.ob.OBAtomAtomIter(a.OBAtom)
            # Check for neighboring atoms in the ring
            n_coords = [pybel.Atom(neigh).coords for neigh in adj if ring.IsMember(neigh)]
            vec1, vec2 = vector(a.coords, n_coords[0]), vector(a.coords, n_coords[1])
            normals.append(np.cross(vec1, vec2))
        # Given all normals of ring atoms and their neighbors, the angle between any has to be 5.0 deg or less
        for n1, n2 in itertools.product(normals, repeat=2):
            arom_angle = vecangle(n1, n2)
            if all([arom_angle > config.AROMATIC_PLANARITY, arom_angle < 180.0 - config.AROMATIC_PLANARITY]):
                return False
        return True

    def find_rings(self):
        """
        Find planar aromatic ring systems
        within molecular object

        Returns
        -------
        rings :class:`list`
            list of namedtuples with all ring
            information extracted i.e.
            atomic coordinates, normal vector
            of ring plane, centroid of ring
            plane, ring size

        ring_centres :class:`list`
            list of tuples with all
            ring centroid coordinates

        References
        ----------
        Definition: aromatic ring definition
        adapted from receptor-ligand interaction
        package as presented in PLIP:
        https://github.com/pharmai/plip
        """
        data = namedtuple('aromatic_ring', 'atoms normal obj center type')
        ring_candidates = self.omol.GetSSSR()

        for ring in ring_candidates:
            ring_atoms = [a for a in self.all_atoms.atoms if ring.IsMember(a.OBAtom)]
            if PiStacks.SMALLEST_RING_SIZE < len(ring_atoms) <= PiStacks.LARGEST_RING_SIZE:
                if ring.IsAromatic() or _is_ring_planar(ring, ring_atoms):
                    ring_type = ring.GetType() if ring.GetType() != '' else 'unknown'

                    ring_atms = [ring_atoms[a].coords for a in [0, 2, 4]]
                    # Probe atoms for normals, assuming planarity
                    ringv1 = vector(ring_atms[0], ring_atms[1])
                    ringv2 = vector(ring_atms[2], ring_atms[0])

                    ring_list = [atom.coords for atom in ring_atoms]

                    self.all_ring_atoms.append(ring_list)

                    self.rings.append(data(atoms=ring_atoms,
                                    normal=unit_vector(np.cross(ringv1, ringv2)),
                                    obj=ring,
                                    center=centroid([ra.coords for ra in ring_atoms]),
                                    type=ring_type))
                    self.ring_centres.append(centroid([ra.coords for ra in ring_atoms]))
        return self.rings, self.ring_centres

    def is_ring_terminal(self, r):
        """ 
        Determines whether the aromatic ring is
        bridged (not-terminal) or terminal.

        Parameters
        ----------
        ring coordinates : :class:`namedtuple`
            The ring to check for its relative
            position in the molecule.


        Returns
        -------
        :class:`str`
            string either saying "terminal"
            or "bridged"

        """
        ring = [atom.coords for atom in r.atoms]
        all_other_rings = [aromatic_ring for aromatic_ring in self.all_ring_atoms if not (aromatic_ring == ring)]
        flat_all_other_rings = list(itertools.chain.from_iterable(all_other_rings))
        common_atoms_counter = sum(1 for x in flat_all_other_rings for y in ring if x == y)
        assert common_atoms_counter == 4 or 2
        return "terminal" if common_atoms_counter == 2 else "bridged"

    def pistacking(self):
        """
        Calculate the pi-pi stacking interactions (parallel and
        orthogonal) for all aromatic rings detected in the
        molecular crystal

        Returns
        -------
        :class:`list`
            list of namedtuples with all
            pi-stacking information

        References
        ----------
        Definition: pi-stacking definition from receptor-
        ligand interaction as presented in PLIP:
        https://github.com/pharmai/plip
        """
        data = namedtuple(
            'pistack', 'ringA ringB distance angle offset \
            type position_ringA position_ringB')

        ringsA = self.rings[:self.rings_per_mol]
        ringsB = self.rings[self.rings_per_mol:]
        for r, l in itertools.product(ringsA, ringsB):
            d = euclidean3d(r.center, l.center)
            b = vecangle(r.normal, l.normal)
            # returns smallest angle between the two ring centres:
            angle = min(b, 180 - b if not 180 - b < 0 else b)

            # calculate slippage (offset) between the two interacting rings
            proj1 = projection(l.normal, l.center, r.center)
            proj2 = projection(r.normal, r.center, l.center)
            offset = min(euclidean3d(proj1, l.center), euclidean3d(proj2, r.center))

            passed = False
            if config.MIN_DIST < d < config.PISTACK_DIST_MAX:
                if 0 < angle < config.PISTACK_ANG_DEV and offset < config.PISTACK_OFFSET_MAX:
                    ptype = 'P' # parallel
                    passed = True
                if (90 - config.PISTACK_ANG_DEV < angle < 90 + config.PISTACK_ANG_DEV) and (offset < config.PISTACK_OFFSET_MAX):
                    ptype = 'T' # tshaped int (orthogonal)
                    passed = True

                if passed:
                    r_pos_ring = self.is_ring_terminal(r)
                    l_pos_ring = self.is_ring_terminal(l)
                    # print("Relevant stack: ", d, "Type of Int: ", ptype, "Offset: ", offset)
                    self.pi_stacks.append(data(ringA=r, ringB=l, distance=d, angle=angle, offset=offset,
                                type=ptype, position_ringA=r_pos_ring,
                                position_ringB=l_pos_ring))
        return self.pi_stacks

    def create_dummy_centroids(self, interacting_mols):
        """
        Create dummy atoms in the position of the centroids
        for each of the aromatic rings.

        !TO DO:
        this method (not) executed, depending on the true/false tag
        as class input argument
        """
        dummy_labelsA = ['H'] * self.rings_per_mol
        dummy_labelsB = ['C'] * self.rings_per_mol
        centroid_labels = dummy_labelsA + dummy_labelsB
        print(centroid_labels, self.rings_per_mol)
        print(self.ring_centres)
        save_files(np.asarray(self.ring_centres), centroid_labels,
                   "centroids", f"{interacting_mols[0]}_{interacting_mols[1]}")