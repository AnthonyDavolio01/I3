import mdtraj as md
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-nc','--numCat',default=6, help='number of ligand CV categories')
parser.add_argument('-s','--start',default=100, help='start frame')
parser.add_argument('-e','--end',default=400, help='end frame')
parser.add_argument('-aa','--numAminos',default=489, help='number of amino acids')
args = parser.parse_args()
aa=int(args.numAminos)
nc=int(args.numCat)
st=int(args.start)
en=int(args.end)
nf=en-st

pos=np.zeros(shape=(aa,nf*nc,3))
for t in range(nc):
    tag=str(t)
    print(tag)
    traj = md.load('output/mc/'+tag+'lig.dcd', top='output/mc/'+tag+'lig.psf')
    for snap in range(st,en):
        print(snap)
        snap2=snap-st
        snap2+=300*t
        for residue in range(0,aa):
            pos[residue][snap2]=traj.xyz[snap][traj.topology.select('resid '+str(residue))].mean(axis=0)

print(pos)
np.save('output/observables/positions.npy',pos)