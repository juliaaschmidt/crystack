# CONSTANTS PER LANDSCAPE:
n_rings=9
#10: 66
#11: 72 
#12: 78
n_atoms=60

# CRYSTAL VARIABLES:
crystal_no=$1

cd res/

# convert .res to .cif
source activate py27
bash convert.sh $crystal_no
conda deactivate

# mv cif to here.
mv opt_"$crystal_no".cif ../

cd ../

# make supercell & compute contacts & append pi-pi stacking degree to results file

conda activate py38contacts
# run execution script with -- var args
python crystack_degree_of_pi_stacking_per_mol_in_crystal.py opt_"$crystal_no".cif $n_atoms $n_rings
# delete dimer_*_*.mol files required for contact analysis
conda deactivate


mkdir dim_"$crystal_no" 
mv real_dimer_*_*.mol dim_"$crystal_no"
zip -r dim_"$crystal_no".zip  dim_"$crystal_no"

rm coms_distance_sorted.csv
#rm real_dimer_*_*.mol
rm -r dim_"$crystal_no"
rm mol_*.xyz
rm nonreal_*.xyz
rm opt_"$crystal_no"_supercell.xyz
rm opt_"$crystal_no".cif
