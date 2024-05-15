#!/bin/bash

# Convert frames to CRD files
#for n in 0 1a 1b 1c 1d 1e 2ab 2ac 2ad 2ae 2bc 2bd 2be 2cd 2ce 2de 3abc 3abd 3abe 3acd 3ace 3ade 3bcd 3bce 3cde 4abcd 4abce 4abde 4acde 4bcde 5 #2ab
#do
#charmm -i code/channel_barrier/charmm/coor$n.inp
#done

# Add sodium to PSF and CRD files
python code/channel_barrier/na_psf.py
python code/channel_barrier/placeNA.py --start 100 --end 4000 -nC 31

# Check profile and narrow location of the barrier
charmm -i code/channel_barrier/charmm/channel_profile_check.inp
python code/channel_barrier/read.py --frames 1 --check_profile --numPoints 180

# Move sodium and calculate energy
for n in 0 1a 1b 1c 1d 1e 2ab 2ac 2ad 2ae 2bc 2bd 2be 2cd 2ce 2de 3abc 3abd 3abe 3acd 3ace 3ade 3bcd 3bce 3cde 4abcd 4abce 4abde 4acde 4bcde 5 #2ab
do
# Alter the start frame in the following. The initial profile was plotted and the script iterates over a smaller range for a speedup.
charmm -i code/channel_barrier/charmm/ion$n.inp
python code/channel_barrier/read.py --tag $n
done
python code/channel_barrier/combine.py -nc 31

# Visually inspect heatmap
python code/channel_barrier/heatmap.py -nc 31