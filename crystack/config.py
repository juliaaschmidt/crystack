
VERBOSE = False  # Set verbose mode
DEBUG = True  # Set debug mode
MAXTHREADS = 1  # Maximum number of main threads for binding site visualization
XML = False
TXT = False
PICS = False
PYMOL = False
STDOUT = False
RAWSTRING = False  # use raw strings for input / output
OUTPATH = './'
BASEPATH = './'
BREAKCOMPOSITE = False  # Break up composite ligands with covalent bonds
ALTLOC = False  # Consider alternate locations
PLUGIN_MODE = False  # Special mode for PLIP in Plugins (e.g. PyMOL)
NOFIX = False  # Turn off fixing of errors in PDB files
NOFIXFILE = False  # Turn off writing to files for fixed PDB structures
PEPTIDES = []  # Definition which chains should be considered as peptide ligands
INTRA = None
KEEPMOD = False
DNARECEPTOR = False
OUTPUTFILENAME = "report"  # Naming for the TXT and XML report files
NOPDBCANMAP = False  # Skip calculation of mapping canonical atom order: PDB atom order

# Configuration file for Protein-Ligand Interaction Profiler (PLIP)
# Set thresholds for detection of interactions
# PISTACK_DIST_MAX = 5.5  # Max. distance for parallel or offset pistacking (McGaughey, 1998)
# PISTACK_ANG_DEV = 30  # Max. Deviation from parallel or perpendicular orientation (in degrees)
# PISTACK_OFFSET_MAX = 2.0  # Maximum offset of the two rings (corresponds to the radius of benzene + 0.5 A)
# PICATION_DIST_MAX = 6.0


# Thresholds for detection (global variables)
BS_DIST = 7.5  # Determines maximum distance to include binding site residues
AROMATIC_PLANARITY = 5.0  # Determines allowed deviation from planarity in aromatic rings
MIN_DIST = 0.5  # Minimum distance for all distance thresholds
# Some distance thresholds were extended (max. 1.0A) if too restrictive too account for low-quality structures
HYDROPH_DIST_MAX = 4.0  # 4.0 # Distance cutoff for detection of hydrophobic contacts
HBOND_DIST_MAX = 4.1  # Max. distance between hydrogen bond donor and acceptor (Hubbard & Haider, 2001) + 0.6 A
HBOND_DON_ANGLE_MIN = 100  # Min. angle at the hydrogen bond donor (Hubbard & Haider, 2001) + 10
PISTACK_DIST_MAX = 5.5  # Max. distance for parallel or offset pistacking (McGaughey, 1998)
PISTACK_ANG_DEV = 30  # Max. Deviation from parallel or perpendicular orientation (in degrees)
PISTACK_OFFSET_MAX = 2.0 # modified to 4.0 (loose) from originally: 2.0 (tight)  # Maximum offset of the two rings (corresponds to the radius of benzene + 0.5 A)
PICATION_DIST_MAX = 6.0  # Max. distance between charged atom and aromatic ring center (Gallivan and Dougherty, 1999)
SALTBRIDGE_DIST_MAX = 5.5  # Max. distance between centers of charge for salt bridges (Barlow and Thornton, 1983) + 1.5
HALOGEN_DIST_MAX = 4.0  # 4.0 Max. distance between oxy. and halogen (Halogen bonds in biological molecules., Auffinger)+0.5
HALOGEN_ACC_ANGLE = 120  # Optimal acceptor angle (Halogen bonds in biological molecules., Auffinger)
HALOGEN_DON_ANGLE = 165  # Optimal donor angle (Halogen bonds in biological molecules., Auffinger)
HALOGEN_ANGLE_DEV = 30  # Max. deviation from optimal angle
WATER_BRIDGE_MINDIST = 2.5  # Min. distance between water oxygen and polar atom (Jiang et al., 2005) -0.1
WATER_BRIDGE_MAXDIST = 4.1  # Max. distance between water oxygen and polar atom (Jiang et al., 2005) +0.5
WATER_BRIDGE_OMEGA_MIN = 71  # Min. angle between acceptor, water oxygen and donor hydrogen (Jiang et al., 2005) - 9
WATER_BRIDGE_OMEGA_MAX = 140  # Max. angle between acceptor, water oxygen and donor hydrogen (Jiang et al., 2005)
WATER_BRIDGE_THETA_MIN = 100  # Min. angle between water oxygen, donor hydrogen and donor atom (Jiang et al., 2005)
METAL_DIST_MAX = 3.0  # Max. distance between metal ion and interacting atom (Harding, 2001)

# Other thresholds
MAX_COMPOSITE_LENGTH = 200  # Filter out ligands with more than 200 fragments
