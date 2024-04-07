from pymol import cmd
import numpy
import argparse
import math
parser = argparse.ArgumentParser()
parser.add_argument('-p','--pdbid',default='6x3z', help='PDB ID / filename prefix')
args=parser.parse_args()
pdbid=args.pdbid


def matriz_inercia(selection):
	'''
	DESCRIPTION

	The method calculates the mass center, the inertia tensor and the eigenvalues and eigenvectors
	for a given selection. Mostly taken from inertia_tensor.py
	'''

	model = cmd.get_model(selection)
	totmass = 0.0
	x,y,z = 0,0,0
	for a in model.atom:
		m = a.get_mass()
		x += a.coord[0]*m
		y += a.coord[1]*m
		z += a.coord[2]*m
		totmass += m
	global cM
	cM = numpy.array([x/totmass, y/totmass, z/totmass])


	I = []
	for index in range(9):
		I.append(0)

	for a in model.atom:
		temp_x, temp_y, temp_z = a.coord[0], a.coord[1], a.coord[2]
		temp_x -= x
		temp_y -= y
		temp_z -= z

		I[0] += a.get_mass() * (temp_y**2 + temp_z**2)
		I[1] -= a.get_mass() * temp_x * temp_y
		I[2] -= a.get_mass() * temp_x * temp_z
		I[3] -= a.get_mass() * temp_x * temp_y
		I[4] += a.get_mass() * (temp_x**2 + temp_z**2)
		I[5] -= a.get_mass() * temp_y * temp_z
		I[6] -= a.get_mass() * temp_x * temp_z
		I[7] -= a.get_mass() * temp_y * temp_z
		I[8] += a.get_mass() * (temp_x**2 + temp_y**2)

	global tensor
	tensor = numpy.array([(I[0:3]), (I[3:6]), (I[6:9])])

	global autoval, autovect, ord_autoval, ord_autovect
	autoval, autovect = numpy.linalg.eig(tensor)
	auto_ord = numpy.argsort(autoval)
	ord_autoval = autoval[auto_ord]
	ord_autovect_complete = autovect[:, auto_ord].T
	ord_autovect = numpy.around(ord_autovect_complete, 3)

	return ord_autovect


cmd.load("input/charmm-gui_input-generator/"+pdbid+"_proa.pdb", "proa")
cmd.load("input/charmm-gui_input-generator/"+pdbid+"_prob.pdb", "prob")
cmd.load("input/charmm-gui_input-generator/"+pdbid+"_proc.pdb", "proc")
cmd.load("input/charmm-gui_input-generator/"+pdbid+"_prod.pdb", "prod")
cmd.load("input/charmm-gui_input-generator/"+pdbid+"_proe.pdb", "proe")
cmd.load("input/charmm-gui_input-generator/"+pdbid+"_heta.pdb", "heta")
cmd.load("input/charmm-gui_input-generator/"+pdbid+"_hetb.pdb", "hetb")



center_of_mass = cmd.centerofmass("all")
cmd.translate([-center_of_mass[0], -center_of_mass[1], -center_of_mass[2]], "all")


z=[0,0,1]
i1=matriz_inercia("all")[0]
ax=numpy.cross(z,i1)
ax = ax / numpy.linalg.norm(ax)
ax = ax.tolist()
ang=numpy.arccos(numpy.dot(z,i1)) * 180 / math.pi
print(i1)
print(ax)
print(ang)
cmd.rotate(ax, angle=ang, origin=[0,0,0], object='all',camera=0)

cmd.save("input/"+pdbid+"_aligned.pdb", "all")
cmd.save("input/"+pdbid+"_aligned_proa.pdb", "proa")
cmd.save("input/"+pdbid+"_aligned_prob.pdb", "prob")
cmd.save("input/"+pdbid+"_aligned_proc.pdb", "proc")
cmd.save("input/"+pdbid+"_aligned_prod.pdb", "prod")
cmd.save("input/"+pdbid+"_aligned_proe.pdb", "proe")

cmd.save("input/"+pdbid+"_aligned_heta.pdb", "heta")
cmd.save("input/"+pdbid+"_aligned_hetb.pdb", "hetb") #hetc in file
exit()
for comp in ['het','pro']:
    for subunit in ['a','b','c','d','e']:
        if comp=='pro' or (comp=='het' and (subunit=='a' or subunit=='b')):
            cmd.load("input/charmm-gui_input-generator/"+pdbid+"_"+comp+subunit+".pdb", pdbid+"_"+comp+subunit)
            cmd.translate([-center_of_mass[0], -center_of_mass[1], -center_of_mass[2]], pdbid+"_"+comp+subunit)
            cmd.rotate(ax, angle=ang, origin=[0,0,0], object='protein',camera=0)
            cmd.save("input/"+pdbid+"_aligned_"+comp+subunit+".pdb", pdbid+"_"+comp+subunit)