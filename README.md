# crystack

A Python module to analyse the intermolecular interactions across an organic
molecular material. This code was written as part of my PhD thesis and the
crystal interactions module.

This repository is of scientific relevance, because it takes a crystallographic
database file structure, which is provided in periodic coordinates, transforms
it to cartesian coordinates and builds a 3x3x3 supercell of the molecular
material. 

Then the most central molecule is selected and all its surrounding molecules
(neighbouring shell) are selected as well, each of them is tested for their
interactions with the central molecule. 

The major focus is on pi-pi stacking interactions which are relevant to the
organic semiconductors industry. Depending on their nature, parallel,
perpendicular etc. they behave differently. Thus, if their distance and angle is
known and their electron transport has been computed elsewhere it can be related
to one another.

## Repository Structure

The crystack folder contains the script to perform the intermolecular analysis.

This repo is structured as follows:

The entire ```crystack``` analysis for a relevant molecular crystal is run via
modifying the automation script: ```run_analysis_for_crystal.sh```

There are several Jupyter notebooks for running the analysis interactively.

The crystal stacking results are stored in a results.csv file and some results
can be found here: ```results/```

Raw Molecular Data: ```data/```

Residuals files of Molecular Crystals: ```res/```

Preliminary results plots: ```Figures/```
 
