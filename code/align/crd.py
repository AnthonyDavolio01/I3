def pdb_to_crd(input_pdb_file, output_crd_file,compsub):
    with open(input_pdb_file, 'r') as pdb_file:
        lines = pdb_file.readlines()

    atom_lines = [line for line in lines if line.startswith('ATOM')]

    with open('input/charmm-gui_input-generator/6dg8_'+comp+subunit+'.crd','r') as orig_file:
        ol = orig_file.readlines()

    for l in range(1,len(ol)):
        ol[l]=ol[l][0:46]+"{:.10f}".format(float(lines[l-1][30:38])).rjust(14)+'      '+"{:.10f}".format(float(lines[l-1][38:46])).rjust(14)+'      '+"{:.10f}".format(float(lines[l-1][46:54])).rjust(14)+ol[l][100::]


    with open(output_crd_file, 'w') as crd_file:
        crd_file.writelines(ol)

for comp in ['het','pro']:
    for subunit in ['a','b','c','d','e']:
        pdb_to_crd('input/6dg8_aligned_'+comp+subunit+'.pdb', 'input/6dg8_aligned_'+comp+subunit+'.crd',comp+subunit)