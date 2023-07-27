"""

AUTOMATED SHORT CONTACTS DIMER CODE

This code reads in the crystal parameters and prints out all the
unique dimers within a chosen *cutoff* (15 A default) away from the ref
molecule in the asymmetric unit. Then for each dimer, it analyses the
intermolecular contacts.

Precautions:

--> Maximum of 15 dimers within NeighShell, but actually you end up with
    less because the equivalent COM-COM distance ones are omitted.

--> Adjust cutoff for different types of molecules (other than helicenes).

--> so far code worked fine on molecules other than helicenes, such as
    pentacene and NTCDI, but use with care on new molecules


Currently employed:

- starting from .cif file and .xyz of asym molecule (without header),
  builds the neighbouring shell --neighbours dimers within a COM-COM
  --cutoff distance.
- All dimers with a duplicate COM-COM distance are omitted.
- tries to detect louvain clusters based on the cosine
  similarity of the dimers to see which dimers are more alike.
  --- EXCHANGE FOR BETTER METHOD ---
- based on the "unique" dimers, calculate all intermolecular short
  contacts possible and print to screen.


Pending actions:

-   find another way of extracting the asymmetric molecule from the
    unit cell than cx1, so that the code can be more widely applicable.
    Avoid the old version of openbabel
-   avoid pywindow
-   compute the energy associated with the contacts? dimer binding energy
    or atomistic picture?

"""
import argparse
import json
from copy import deepcopy

from interactions import Interactions

Sunshine = "#f5cd7e"
Coral = "#eb6b67"
Silver = "#e1dcdc"
Peacock = "#6b9ba5"


    # avDist is the average distance of the terminal rings midpoint to the
    # dimer centre-of-mass (int vs b2b). DistDiff is the difference in
    # the midpoint of terminal rings distance to the dimer centre-of-mass
    # (trans vs. herringbone) mols_orient is the orientation angle of the
    # molecules, -1 perfectly antiparallel, 0 orthogonal and
    # 1 perfectly parallel.
    # >0.95 translational/herringbone
    # < 0 back2back
    # < -0.4 interlocked
    # NOTES TO MYSELF:
    # interlocked features; many hydrophobic interactions, no pi stacks,
    # sometimes one short contact
# ----------------------------------------------------------------------------#


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description= "Generate all unique \
                                     Dimer files of neighbouring shell \
                                     within a COM-COM cutoff. ")
    parser.add_argument('filns', type=str, help='List of all filenames \
                        to be compared')
    # parser.add_argument('cell_lengths', metavar='cell_lengths', \
    #                     type=float, nargs=3, help='cell lengths a,b,c')
    # parser.add_argument('cell_angles',metavar='cell_angles', \
    #                     type=float, nargs=3, help='alpha, beta, gamma')
    # parser.add_argument('sg', type=str, help='Crystal Spacegroup')
    # parser.add_argument('-cutoff','--cutoff', default=15.0,
    #                     help='Set individual cutoff (default: 15.0 A)')
    # parser.add_argument('-neighbours', '--neighbours',type=int,
    #                     default=15, help='Set individual number of \
    #                     nearest neighbours (default: 15)')
    parser.add_argument('-vdw-plus-extra', '--vdwextra', type=float,
                        default=0., help='Set the intermolecular short \
                        contacts threshold, default is sum of the vdw \
                        radii plus +2.0 A (Mercury Default +0.0 A )')
    # parser.add_argument('outn',type=str,
    #                     help='Output name - crystal structure identifier)')
    args = parser.parse_args()

    filns = [line.rstrip('\n') for line in open(args.filns)]
    # _________________STUDY INTERACTIONS IN EACH DIMER _______________________

    motifs_dict = {}

    for i, filn in enumerate(filns):
        # for key, value in df.iteritems():
        """ Show pi-pi stacks within dimer """
        print("________________________________________________")
        print(filn)
        #red_name = filn.split("real_")[-1]

        dimer = Interactions(f"{filn}")
        dimer.make_mol_obj()

        #dimer.calculate_pi_stacks()  # 
        dimer.calculate_all_contacts(vdwextra=args.vdwextra)

        print("________________________________________________")
    
        pi_stack_info = [(item.distance, item.type, item.offset, item.position_ringA, item.position_ringB) for item in dimer.pi_stacks]
        short_contacts_info = [(item.typeA, item.typeB, item.distance) for
        item in dimer.short_stacks]
        hydrophobic_ints_info = [(item.typeA, item.typeB, item.distance) for item in dimer.hydro_stacks]
        # pi_stack_d = [item.distance for item in dimer.pi_stacks]
        # pi_stack_offset = [item.offset for item in dimer.pi_stacks] 
        #dimer.create_dummy_centroids()

        print(pi_stack_info)
        # print(short_contacts_info)
        print("________________________________________________")

        dimer.__repr__()
        with open(f'{red_name}.json', 'w') as f:
            pass
            # pi_stack = {}
            # pi_stack["pi_stack"] = pi_stack_d
            # short_contacts = {}
            # short_contacts["short_contacts"] = short_contacts_d
            # hydrophobic_ints = {}
            # hydrophobic_ints["hydrophobic_interactions"] = hydrophob_d
            # orientation = {}
            # orientation['orientation'] = mols_orient
            # json.dump(pi_stack, f)
            # json.dump(short_contacts, f)
            # json.dump(hydrophobic_ints, f)
            # json.dump(mols_orient)

        # visualise centroids?
    # analyse_contacts("Dimer_1", fake_centroids_xyz=False)
