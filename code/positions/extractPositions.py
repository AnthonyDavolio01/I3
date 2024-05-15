import mdtraj as md
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-nc','--numCat',default=6, help='number of ligand CV categories')
parser.add_argument('-s','--start',default=100, help='start frame')
parser.add_argument('-e','--end',default=4000, help='end frame')
parser.add_argument('-aa','--numAminos',default=489*5, help='number of amino acids')
args = parser.parse_args()
aa=int(args.numAminos)
nc=int(args.numCat)
st=int(args.start)
en=int(args.end)
nf=en-st

#pos=np.zeros(shape=(aa,nf*nc,3))
pos=np.zeros(shape=(aa,nf*nc,3))
#pos=np.load('output/observables/positions.npy')


for t in range(nc):
    tag=str(t)
    print('category '+tag)
    for snap in range(st,en):
        print('snapshot '+str(snap))
        ind=0
        with open('output/observables/barrier/coors/'+tag+'_lig.'+str(snap)+'.coor','r') as file:
            coors=file.readlines()
            a=[]
            aNum=1
            for line in coors[4::]:
                l=line.split()
                if l[7]=='PROA' or l[7]=='PROB' or l[7]=='PROC' or l[7]=='PROD' or l[7]=='PROE':
                    if int(l[1])==aNum:
                        a.append([float(l[4]),float(l[5]),float(l[6])])
                    else:
                        p=np.asarray(a).mean(axis=0)
                        pos[ind][t*(en-st)+snap-st]=p
                        a=[]
                        ind+=1
                        aNum=int(l[1])
            p=np.asarray(a).mean(axis=0)
            pos[ind][t*(en-st)+snap-st]=p

                
print(pos)



np.save('output/observables/positions.npy',pos)
exit()

# Slower, using mdtraj
for t in range(2,nc):#nc
    tag=str(t)
    print('category '+tag)
    traj = md.load('output/mc/'+tag+'lig.dcd', top='output/mc/'+tag+'lig.psf')
    for snap in range(st,en):
        print('snapshot '+str(snap))
        ind=0
        for chain_idx, chain in enumerate(traj.topology.chains):
            for residue in chain.residues:
                p=traj.xyz[snap][traj.topology.select('resid '+str(residue.index)+' and chainid '+str(chain_idx))].mean(axis=0)
                print(p)
                pos[ind][t*(en-st)+snap-st]=p
                ind+=1
                
print(pos)




np.save('output/observables/positions.npy',pos)
exit()
#for snap in range(st,en):
#        print(snap)
#        snap2=snap-st
#        snap2+=(en-st)*t
#        for residue in range(0,aa):
#            pos[residue][snap2]=traj.xyz[snap][traj.topology.select('resid '+str(residue)+' and chainid '+str(cid))].mean(axis=0)

print(pos)
np.save('output/observables/positions.npy',pos)