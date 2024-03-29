import numpy as np
from sklearn.decomposition import PCA
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-nc','--numCat',default=6, help='number of ligand CV categories')
parser.add_argument('-f','--frames',default=300, help='number of frames after equilibration')
parser.add_argument('-aa','--numAminos',default=489, help='number of amino acids')
args = parser.parse_args()
f=int(args.frames)
nc=int(args.numCat)
naa=int(args.numAminos)

positions=np.load('output/observables/positions.npy')
aCV=np.load('output/observables/barrier/channel/all_barrier.npy')
l=f*nc

ligCV = []
for c in range(nc):
    for fr in range(f):
        ligCV.append(c)
ligCV = np.array(ligCV)

bCV=ligCV

numBins=10

Abins=[]
for i in range(numBins):
    Abins.append((np.max(aCV)-np.min(aCV))*i/numBins+np.min(aCV))

aD=np.digitize(aCV, bins=Abins)-1
bD=bCV.astype(int)




#PCA on positions, use first eigenvector
positions_firstEig=np.zeros((naa,l))
for aa in range(naa):
    pca = PCA(n_components=1)
    pca.fit(positions[aa])
    positions_firstEig[aa]=pca.transform(positions[aa]).T[0].T


zD=np.zeros((naa,l))
for aa in range(naa):
    paa=positions_firstEig[aa]
    zBins=[]
    for i in range(numBins):
        zBins.append((np.max(paa)-np.min(paa))*i/numBins+np.min(paa))
    zDaa=np.digitize(paa, bins=zBins)-1
    zD[aa]=zDaa

zD=zD.astype(int)


#abDist=np.zeros((numBins,6))
abzDist=np.zeros((numBins,nc,naa,numBins))
for i in range(l):
    for aa in range(naa):
        abzDist[aD[i]][bD[i]][aa][zD[aa][i]]+=1/l

i2xy=np.zeros(naa)
i2xyz=np.zeros(naa)

for x,y,aa,z in np.ndindex(abzDist.shape):
    if abzDist[x,y,aa,z]!=0:
        i2xy[aa]+=abzDist[x,y,aa,z]*np.log(np.sum(abzDist[x,y,aa,:])/(np.sum(abzDist[x,:,aa,:])*np.sum(abzDist[:,y,aa,:])))
        i2xyz[aa]+=abzDist[x,y,aa,z]*np.log((np.sum(abzDist[x,y,aa,z])*np.sum(abzDist[:,:,aa,z]))/(np.sum(abzDist[x,:,aa,z])*np.sum(abzDist[:,y,aa,z])))

i3=i2xy-i2xyz

np.save('output/I3.npy',i3)

print(i3)


min_val = i3.min()
max_val = i3.max()
normalized_arr = (i3 - min_val) / (max_val - min_val)


numColor=50
sorted_indices = np.argsort(normalized_arr)[::-1]

# Get the indices of the 50 highest values
top_50_indices = sorted_indices[:numColor]


with open("output/i3colors.txt", "w") as file:
    for aa in top_50_indices:
        file.write("set_color me"+str(aa)+", [0,"+str(1-normalized_arr[aa])+","+str(normalized_arr[aa])+"]\n")

with open("output/i3set.txt", "w") as file:
    for aa in top_50_indices:
#        file.write("set_color "+str(aa)+", [0,0,"+str(normalized_arr[aa])+"]\n")
        file.write("cmd.color(\'me"+str(aa)+"\', \'resi "+str(aa+8)+"\')\n")