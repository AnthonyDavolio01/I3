import numpy as np
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('-f','--frames',default=300, help='number of frames after equilibration')
parser.add_argument('-nc','--numCat',default=6, help='number of ligand CV categories')

args = parser.parse_args()

nf=int(args.frames)
nc=int(args.numCat)
dim=nf*nc
barrier=np.zeros(dim)

for t in range(nc):
    tag=str(t)
    temp=np.load('output/observables/barrier/channel/barrier'+tag+'.npy')
    barrier[nf*t:nf*(t+1)]=temp

np.save('output/observables/barrier/channel/all_barrier.npy',barrier)
print(barrier.shape)
