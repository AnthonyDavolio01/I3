
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('-s', '--start', default=100, help='The equilibrated frame')
parser.add_argument('-e', '--end', default=4000, help='The end frame')
parser.add_argument('-nC', '--numCat', default=31, help='number of categories')
    
args = parser.parse_args()
cat=['0', '1a', '1b', '1c', '1d', '1e', '2ab', '2ac', '2ad', '2ae', '2bc', '2bd', '2be', '2cd', '2ce', '2de', '3abc', '3abd', '3abe', '3acd', '3ace', '3ade', '3bcd', '3bce', '3cde', '4abcd', '4abce', '4abde', '4acde', '4bcde', '5']


for t in cat:
    tag=str(t)
    print(tag)
    for i in range(int(args.start),int(args.end)):
        with open('output/observables/barrier/coors/'+tag+'_5ht3.'+str(i)+'.coor', 'r') as f:
            lines = f.readlines()
            atomInd=int(lines[-1].split()[0])+1
            resInd=int(lines[-1].split()[1])+1
            lines[3] = '     '+str(atomInd)+'  EXT\n'
            lines.append('     '+str(atomInd)+'      '+str(resInd)+'  SOD       SOD            00.0000000000       00.0000000000      -110.0000000000  HETF      6               0.0000000000')

        with open('output/observables/barrier/coors/'+tag+'_na_5ht3.'+str(i)+'.coor', 'w') as g:
            g.writelines(lines)