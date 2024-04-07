#!/bin/bash

# Convert frames to CRD files
for n in 0 1 2
do
charmm -i code/channel_barrier/charmm/coor$n.inp
done

# Add ion (chlorine) to PSF and CRD files
python code/channel_barrier/na_psf.py
python code/channel_barrier/placeNA.py --start 100 --end 400 -nC 3

# Check profile and narrow location of the barrier
charmm -i code/channel_barrier/charmm/channel_profile_check.inp
python code/channel_barrier/read.py --frames 1 --check_profile --numPoints 100

# Move ion and calculate energy
for n in 0 1 2
do
# Alter the start frame in the following. The initial profile was plotted and the script iterates over a smaller range for a speedup.
charmm -i code/channel_barrier/charmm/ion$n.inp
python code/channel_barrier/read.py --tag $n
done
python code/channel_barrier/combine.py -nc 3

# Visually inspect heatmap
python code/channel_barrier/heatmap.py -nc 3
# Barrier raises on one ligand bc assymetrical.