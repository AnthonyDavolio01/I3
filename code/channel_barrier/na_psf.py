
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('-nC', '--numCat', default=3, help='number of categories')
    
args = parser.parse_args()

for n in range(int(args.numCat)):
    with open(f'output/mc/{n}lig.psf','r') as start_psf:
        l1=start_psf.readlines()
        natom=int(l1[6].split()[0])+1
        new6=str(natom).rjust(10)+' !NATOM\n'
        l1[6]=new6
        l1[natom+6]=str(natom).rjust(10)+" HETF     6        CLA      CLA      CLA       1.00000       22.9898           0\n\n"

    with open(f'output/observables/barrier/psf/{n}lig.psf','w') as na_psf:
        na_psf.writelines(l1)




#"     38781 HETF     6        SOD      SOD      SOD       1.00000       22.9898           0"