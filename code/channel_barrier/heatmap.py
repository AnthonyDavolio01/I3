import matplotlib.pyplot as plt
import numpy as np
import argparse

data=np.load('output/observables/barrier/channel/all_barrier.npy')
parser = argparse.ArgumentParser()
parser.add_argument('-nc','--numCat',default=3, help='number of ligand CV categories')
parser.add_argument('-f','--numFrames',default=300, help='number of frames after equilibration')
args = parser.parse_args()

ligCV = []
for c in range(int(args.numCat)):
    for fr in range(int(args.numFrames)):
        ligCV.append(c)
ligCV = np.array(ligCV)

heatmap, xedges, yedges = np.histogram2d(ligCV, data, bins=20)

xmesh, ymesh = np.meshgrid(xedges[:-1], yedges[:-1])

plt.pcolormesh(xmesh, ymesh, heatmap.T, cmap='viridis')
plt.colorbar(label='Density')

plt.xlabel('# of Ligands Bound')
plt.ylabel('Channel Energy Barrier')
plt.title('Energy Barrier vs # of Ligands Bound Heat Density Map')

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