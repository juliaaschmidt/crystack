"""
This input-output (IO) module contains all functionality required to
read in files and generate output files.

* TO DO:
- dummy atoms supported
- real/nonreal xyz file flags
- generated real_xyz file headers

"""
import numpy as np
import subprocess as sp
import pandas as pd

# class Input(object):

#     def __init__(self):
#         self._load_funcs = {
#             '.xyz': self._read_xyz,
#             '.mol': self._read_mol,
#         }




def save_files(file_to_save, labels, filn, mol_number):
    """ Saves all the dimers created to a separate xyz file.
    ! TO DO: HEADER FILE FOR REAL-xyz generation 
    """
    trans = np.array(list(zip(labels, file_to_save[:,0], file_to_save[:,1],\
        file_to_save[:,2])), dtype=[('labels','U8'),('trans[:,0]',float),\
        ('trans[:,1]', float), ('trans[:,2]', float)])

    out_filn = f"{filn}_{mol_number}.xyz"
    moln = f"{filn}_{mol_number}.mol"

    # with open(out_filn, 'a') as f:
    #     f.write(str(len(labels)) + "\n\n")
        
    np.savetxt(f"real_{out_filn}", trans, delimiter=" 	",\
    fmt=["%s"]+["%f"]+["%f"]+["%f"])


    # print(f"Saving {out_filn} to file.")
    # sp.call(f"cat header {out_filn} > real_{out_filn}", shell=True)
    sp.call(f"obabel -ixyz real_{out_filn} -omol > {moln}", shell=True)


def get_mol_centroids(self, stringlist):
    """
    Function to get the centroids of the two interacting molecules
    A and B based on the coms_distance_sorted.csv file.
    """
    centroids_list = pd.read_csv("coms_distance_sorted.csv", skiprows=1, dtype=float, names=['id','x','y','z','dist'])
    centroids_list['id'] = centroids_list["id"].astype('int32')
    # print(stringlist)
    # print(centroids_list)
    mola = centroids_list[centroids_list["id"] == int(stringlist[0])]
    molb = centroids_list[centroids_list["id"] == int(stringlist[1])]
    centreA_col = mola[['x','y','z']]# df
    centreB_col = molb[['x','y','z']]

    centreA = centreA_col.iloc[0].values #row of df
    centreB = centreB_col.iloc[0].values
    print("Centroid position of molecules in dimer: ")
    print("Centre A is: ", centreA)
    print("Centre B is: ", centreB)
    return centreA, centreB

