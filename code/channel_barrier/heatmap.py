import matplotlib.pyplot as plt
import numpy as np
import argparse

data=np.load('output/observables/barrier/channel/all_barrier.npy')
parser = argparse.ArgumentParser()
parser.add_argument('-nc','--numCat',default=31, help='number of ligand CV categories')
parser.add_argument('-f','--numFrames',default=3900, help='number of frames after equilibration')
args = parser.parse_args()

ligCV = []
for c in range(int(args.numCat)):
    for fr in range(int(args.numFrames)):
        ligCV.append(c)
ligCV = np.array(ligCV)

l2=[]
d2=[]
for i in range(len(data)):
    if data[i]!=0:
        l2.append(ligCV[i])
        d2.append(data[i])

print(len(ligCV))
print(len(l2))


heatmap, xedges, yedges = np.histogram2d(l2, d2, bins=20)

xmesh, ymesh = np.meshgrid(xedges[:-1], yedges[:-1])

plt.pcolormesh(xmesh, ymesh, heatmap.T, cmap='viridis')
plt.colorbar(label='# of snapshots')

plt.xlabel('# of Ligands Bound')
plt.ylabel('Channel Energy Barrier (KCal / Mol)')
plt.title('Energy Barrier against # of Ligands Bound Heat Density Map')

plt.show()

exit()

data = data.reshape(int(args.numCat), -1)
for lcv in data:
    count=0
    for f in lcv:
        if f<29:
            count+=1
    print(count)
exit()