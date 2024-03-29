# Allosteric Site Searching with Conditional Mutual Information

### Requirements:
 - CHARMM c48b2
 - PyMol 2.5+
 - matplotlib
 - mdtraj==1.9.9
 - numpy
 - sklearn


## Step 1: Generate Inputs
Use CHARMM-GUI PDB reader & manipulator to generate the input files.
https://www.charmm-gui.org/?doc=input/pdbreader

We used 6DG8 from RCSB and translated its center of mass to the origin with the following script using PyMol.
```bash
python code\align\align.py
python code\align\crd.py
```

## Step 2: Run CHARMM Monte Carlo
Generate the monte carlo ensembles for all six ligand categories:
```bash
bash code/monte_carlo/run_MC.sh
```
This will run the six monte carlo simulations in parallel.

## Step 3: Obtain Observables
Check the energy autocorrelation function:
```bash
charmm -i code/energy_acf/acf.inp
python code/energy_acf/acfPlot.py
```
the equilibrated frame is the start frame for the next commands.

Calculate the Ion Channel Barrier of each frame:
```bash
bash code/channel_barrier/calc_barrier.sh
```

Extract the position of amino acids:
```bash
python extractPositions.py
```

## Step 4: 
Calculate Conditional Mutual Information:
```bash
python code/I3.py
```

### Visualization
Open the input/6dg8_aligned.pdb with PyMol
Paste the contents of output/i3colors.txt and output/i3set.txt directly into PyMol