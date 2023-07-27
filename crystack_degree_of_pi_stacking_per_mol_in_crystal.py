#!/usr/bin/env python
# coding: utf-8

# In[1]:


from crystack.interactions import Interactions
from crystack import cif2supercell
import re
import itertools
import numpy as np
import sys


# In[2]:


# Crystal specifics
cif_name = sys.argv[1] #"opt_1.cif"
n_atoms = int(sys.argv[2]) #78
n_rings = int(sys.argv[3]) #12

#print(n_atoms, n_rings)
# CONSTANTS
neighshell_size = 15

crystal_no = re.findall("[0-9]+",cif_name)[0]
orient_dict = cif2supercell.run_dimer_generator(cif_name, n_atoms, out_mol=neighshell_size)


# In[3]:


from glob import glob
dimer_files = glob('real_dimer*mol')


# In[4]:


import re 
pi_stacks = []
pi_stack_dir = {}
dominant_dir = {}
for f in dimer_files: 
    basen = f.split(".")[0]
    last_val = re.findall("[0-9]+", basen)
    mol_int = int(last_val[-1])
    dominant_dir[mol_int] = orient_dict[mol_int]
    dimer = Interactions(f"{basen}")
    dimer.make_mol_obj()
    pi_stack_dimer = dimer.calculate_pi_stacks()
    pi_stacks.append(pi_stack_dimer)
    pi_stack_dir[mol_int] = pi_stack_dimer 


# In[5]:


pi_dists = []
pi_angles = []
pi_types = []
dominant_dir_val = []
interacting_mols = []
unique_rings_on_mol = []
for i, dimer_no in enumerate(pi_stack_dir):
    for i, n in enumerate(pi_stack_dir[dimer_no]):
        if (n.type == "P"):
            #print(f"{n.type} only.")
            unique_rings_on_mol.append(n.ringA.center)
            pi_dists.append(n.distance)
            pi_angles.append(n.angle)
            #pi_types.append(n.type)
            dominant_dir_val.append(dominant_dir[dimer_no])
            interacting_mols.append(dimer_no)


# In[6]:


crystal_no = re.findall("[0-9]+",cif_name)[0]


# In[7]:


import itertools
import numpy as np

m = unique_rings_on_mol
n = [list(i) for i in set(map(tuple, m))]
n_rings_involved = len(n)
frac = n_rings_involved/n_rings
print(f"Final score is {n_rings_involved} out of {n_rings} rings are involved ({np.round(frac,2)}%).")


# In[8]:


print(f"{crystal_no},{n_rings_involved},{n_rings},{np.round(frac,2)}")


# In[9]:


with open("results.csv", "a") as f:
    f.write(f"{crystal_no},{n_rings_involved},{n_rings},{np.round(frac,2)}\n")


# In[ ]:




