import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-p','--pdbid',default='6x3z', help='PDB ID / filename prefix')
args=parser.parse_args()
pdbid=args.pdbid

def pdb_to_crd(input_pdb_file, output_crd_file,compsub):
    with open(input_pdb_file, 'r') as pdb_file:
        lines = pdb_file.readlines()

    atom_lines = [line for line in lines if line.startswith('ATOM')]

    with open('input/charmm-gui_input-generator/'+pdbid+'_'+comp+subunit+'.crd','r') as orig_file:
        ol = orig_file.readlines()

    for l in range(1,len(ol)):
        ol[l]=ol[l][0:46]+"{:.10f}".format(float(lines[l-1][30:38])).rjust(14)+'      '+"{:.10f}".format(float(lines[l-1][38:46])).rjust(14)+'      '+"{:.10f}".format(float(lines[l-1][46:54])).rjust(14)+ol[l][100::]


    with open(output_crd_file, 'w') as crd_file:
        crd_file.writelines(ol)

for comp in ['het','pro']:
    for subunit in ['a','b','c','d','e']:
        if comp=='pro' or (comp=='het' and (subunit=='a' or subunit=='b')):
            pdb_to_crd('input/'+pdbid+'_aligned_'+comp+subunit+'.pdb', 'input/'+pdbid+'_aligned_'+comp+subunit+'.crd',comp+subunit)