import numpy as np
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('-f','--frames',default=3900, help='number of frames after equilibration')
parser.add_argument('-nc','--numCat',default=31, help='number of ligand CV categories')

args = parser.parse_args()

nf=int(args.frames)
nc=int(args.numCat)
dim=nf*nc
barrier=np.zeros(dim)
cat=['0', '1a', '1b', '1c', '1d', '1e', '2ab', '2ac', '2ad', '2ae', '2bc', '2bd', '2be', '2cd', '2ce', '2de', '3abc', '3abd', '3abe', '3acd', '3ace', '3ade', '3bcd', '3bce', '3cde', '4abcd', '4abce', '4abde', '4acde', '4bcde', '5']

spot=0
for t in cat:
    tag=str(t)
    temp=np.load('output/observables/barrier/channel/barrier'+tag+'.npy')
    barrier[nf*spot:nf*(spot+1)]=temp
    spot+=1

np.save('output/observables/barrier/channel/all_barrier.npy',barrier)
print(barrier.shape)
