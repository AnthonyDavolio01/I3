
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('-nC', '--numCat', default=31, help='number of categories')

args = parser.parse_args()

cat=['0', '1a', '1b', '1c', '1d', '1e', '2ab', '2ac', '2ad', '2ae', '2bc', '2bd', '2be', '2cd', '2ce', '2de', '3abc', '3abd', '3abe', '3acd', '3ace', '3ade', '3bcd', '3bce', '3cde', '4abcd', '4abce', '4abde', '4acde', '4bcde', '5']


for n in cat:
    with open(f'output/mc/{n}lig.psf','r') as start_psf:
        l1=start_psf.readlines()
        natom=int(l1[6].split()[0])+1
        new6=str(natom).rjust(10)+' !NATOM\n'
        l1[6]=new6
        l1[natom+6]=str(natom).rjust(10)+" HETF     6        SOD      SOD      SOD       1.00000       22.9898           0\n\n"

    with open(f'output/observables/barrier/psf/{n}lig.psf','w') as na_psf:
        na_psf.writelines(l1)




#"     38781 HETF     6        SOD      SOD      SOD       1.00000       22.9898           0"