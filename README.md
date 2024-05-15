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

Extract the energy of amino acids:
```bash
for n in 0 1a 1b 1c 1d 1e 2ab 2ac 2ad 2ae 2bc 2bd 2be 2cd 2ce 2de 3abc 3abd 3abe 3acd 3ace 3ade 3bcd 3bce 3cde 4abcd 4abce 4abde 4acde 4bcde 5
do
charmm -i code/nrgs/nrgs$n.inp
done
python code/nrgs/read.py ### at this step, still need to fix
```


## Step 4: Calculate Thermodynamic Coupling Function
```bash
```

### -- Visualization --
Open the input/6dg8_aligned.pdb with PyMol,
paste the contents of output/i3colors.txt then output/i3set.txt directly into PyMol.


### -- Deprecated --
Extract the position of amino acids:
```bash
python code/positions/extractPositions.py
``` 
Calculate Conditional Mutual Information
```bash
python code/I3.py
```


L:
cd repositories
cd I3
wsl -u root