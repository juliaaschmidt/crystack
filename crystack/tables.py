#!/usr/bin/env python3
"""
The table module lists general-purpose chemical data,
required for the interaction analysis

Sources:
    1. Data taken from stk.

Note:
- vdw radii
- atomic number
"""


atom_vdw_radii = {
        'Al': 2, 'Sb': 2, 'Ar': 1.88, 'As': 1.85, 'Ba': 2,
        'Be': 2, 'Bi': 2, 'B': 2, 'Br': 1.85, 'Cd': 1.58,
        'Cs': 2, 'Ca': 2, 'C': 1.7, 'Ce': 2, 'Cl': 1.75,
        'Cr': 2, 'Co': 2, 'Cu': 1.4, 'Dy': 2, 'Er': 2,
        'Eu': 2, 'F':  1.47, 'Gd': 2, 'Ga': 1.87, 'Ge': 2,
        'Au': 1.66, 'Hf': 2, 'He': 1.4, 'Ho': 2, 'H': 1.09,
        'In': 1.93, 'I': 1.98, 'Ir': 2, 'Fe': 2, 'Kr': 2.02,
        'La': 2, 'Pb': 2.02, 'Li': 1.82, 'Lu': 2, 'Mg': 1.73,
        'Mn': 2, 'Hg': 1.55, 'Mo': 2, 'Nd': 2, 'Ne': 1.54,
        'Ni': 1.63, 'Nb': 2, 'N':  1.55, 'Npl':  1.55, 'Os': 2,
        'O': 1.52,
        'Pd': 1.63, 'P': 1.8, 'Pt': 1.72, 'K': 2.75, 'Pr': 2,
        'Pa': 2, 'Re': 2, 'Rh': 2, 'Rb': 2, 'Ru': 2, 'Sm': 2,
        'Sc': 2, 'Se': 1.9, 'Si': 2.1, 'Ag': 1.72, 'Na': 2.27,
        'Sr': 2, 'S': 1.8, 'Ta': 2, 'Te': 2.06, 'Tb': 2,
        'Tl': 1.96, 'Th': 2, 'Tm': 2, 'Sn': 2.17, 'Ti': 2,
        'W': 2, 'U':  1.86, 'V':  2, 'Xe': 2.16, 'Yb': 2,
        'Y': 2, 'Zn': 1.29, 'Zr': 2, 'X':  1.0, 'D':  1.0,
        'O2': 1.52,
}

atomic_num = {
        "H": 1, "He": 2, "Li":3, "Beryllium":4, "Boron":5,
        "C": 6, "N": 7, "O":8, "F":9, "Ne":10, "Na":11,
        "Mg": 12, "Al":13, "Si":14, "P":15, "S":16,
        "Cl": 17, "Ar":18, "K":19, "Ca":20, "Sc":21, "Ti":22,
        "V": 23, "Cr": 24, "Mn":25, "Fe": 26, "Co":27, "Ni":28,
        "Cu":29, "Zn":30, "Ga":31, "Ge": 32, "As":33, "Br":35,
        "Se":34, "Kr":36, "Rb":37, "Sr": 38, "Y":39, "Zr":40,
        "Nb":41, "Mo":42, "Tc":43, "Ru": 44, "Rh":45, "Pd":46,
        "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "I": 53
}
