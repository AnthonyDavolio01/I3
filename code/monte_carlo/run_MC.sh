#!/bin/bash

### Serial

# # # Sometimes the monte carlo will unfold the protein.
# # # In this case, rerun that simulation under a different MC ISEED

## Loop from n=0 to n=5
#for n in 0 1 2 3 4 5
#do
#    charmm -i code/monte_carlo/lig$n.inp
#done


### Parallel

# Define a function for running charmm

run_charmm() {
    charmm -i "$1"
}
#
## Export the function to make it available to GNU Parallel
export -f run_charmm
#
## Run commands in parallel using GNU Parallel
seq 0 5 | parallel -j "$(nproc)" run_charmm code/monte_carlo/lig{}.inp


# NOTE: 
# If you are modifying this code to probe your own protein, 
# remember to alter the disulfide bonds in the 
# monte_carlo lig.inp scripts and read in your own coordinate files.
