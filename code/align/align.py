from pymol import cmd



# The long axis was already aligned to the z-axis in the experimental coordinates.


cmd.load("input/charmm-gui_input-generator/step1_pdbreader.pdb", "protein")
center_of_mass = cmd.centerofmass()

cmd.translate([-center_of_mass[0], -center_of_mass[1], -center_of_mass[2]], "protein")
cmd.save("input/6dg8_aligned.pdb", "protein")




for comp in ['het','pro']:
    for subunit in ['a','b','c','d','e']:
        cmd.load("input/charmm-gui_input-generator/6dg8_"+comp+subunit+".pdb", "6dg8_"+comp+subunit)
        cmd.translate([-center_of_mass[0], -center_of_mass[1], -center_of_mass[2]], "6dg8_"+comp+subunit)
        cmd.save("input/6dg8_aligned_"+comp+subunit+".pdb", "6dg8_"+comp+subunit)