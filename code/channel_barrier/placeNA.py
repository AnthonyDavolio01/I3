
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('-s', '--start', default=100, help='The equilibrated frame')
parser.add_argument('-e', '--end', default=400, help='The end frame')
parser.add_argument('-nC', '--numCat', default=3, help='number of categories')

args = parser.parse_args()


for t in range(0,int(args.numCat)):
    tag=str(t)
    print(tag)
    for i in range(int(args.start),int(args.end)):
        with open('output/observables/barrier/coors/'+tag+'_lig.'+str(i)+'.coor', 'r') as f:
            lines = f.readlines()
            atomInd=int(lines[-1].split()[0])+1
            resInd=int(lines[-1].split()[1])+1
            lines[3] = '     '+str(atomInd)+'  EXT\n'
            lines.append('     '+str(atomInd)+'      '+str(resInd)+'  CLA       CLA             7.0000000000       -7.0000000000      -50.0000000000  HETF      6               0.0000000000')

        with open('output/observables/barrier/coors/'+tag+'_na_lig.'+str(i)+'.coor', 'w') as g:
            g.writelines(lines)