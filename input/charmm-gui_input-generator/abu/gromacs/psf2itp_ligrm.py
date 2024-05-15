"""
Generated by CHARMM-GUI (http://www.charmm-gui.org)

psf2itp.py

This program is for GROMACS input generation with CHARMM36 FF.

Correspondance: jul316@lehigh.edu or wonpil@lehigh.edu
Last update: August 20, 2015
"""

import sys
import os

# conversion constant between kcal and kJ
kcal2kJ = 4.184

# specification of character used for comments in charmm ff file
comment = '!'

if len(sys.argv) != 3:
    print "useage: psf2itp.py topparPath psfFile"
    exit()

#-------------------------------------------#
#               Read parFiles               #
#-------------------------------------------#

path = sys.argv[1]

parFiles = []
for filename in os.listdir(path):
    fileext = filename.split('.')[-1]
    if fileext == 'prm' or fileext == 'rtf' or fileext == 'str':
        parFiles.append(path+'/'+filename)

if len(parFiles) == 0:
    print "Fatal error: incorrect toppar path"
    exit()

index2type     = {}
type2mass      = {}
type2bond      = {}
type2angle     = {}
type2dihedral  = {}
type2mult      = {}
type2improper  = {}
type2nonbonded = {}
type2pair      = {}
type2nbfix     = {}

dihWilds = []
impWilds = []

cmap      = []
cmapCount = []
cmapData  = []

for parfile in parFiles:
    readprm = None
    line_continue = None
    for line in open(parfile, 'r'):
        if line.find(comment) >= 0: line = line.split(comment)[0]
        line = line.upper().strip()
        if len(line) > 0:
            if line_continue:
                line = line_continue + line
                line_continue = None
            if line.endswith('-'):
                line_continue = line[:-1]
                continue
            if line.startswith('RESI'):
                readprm = None
                continue
            if line.startswith('PRES'):
                readprm = None
                continue
            if line.startswith('ATOMS'):
                readprm = 'atom'
                continue
            if line.startswith('BONDS'):
                readprm = 'bond'
                continue
            if line.startswith('ANGLES'):
                readprm = 'angle'
                continue
            if line.startswith('DIHEDRALS'):
                readprm = 'dihedral'
                continue
            if line.startswith('IMPROPER'):
                readprm = 'improper'
                continue
            if line.startswith('CMAP') and line.endswith('CMAP'):
                readprm = 'cmap'
                continue
            if line.startswith('NONBONDED'):
                readprm = 'nonbonded'
                continue
            if line.startswith('NBFIX'):
                readprm = 'nbfix'
                continue
            if line.startswith('HBOND'):
                readprm = None
                continue
            if line.upper().startswith('END'):
                readprm = None
                continue

            segments = line.split()

            if readprm == 'atom':
                type2mass[segments[2]]  = float(segments[3])
                index2type[segments[1]] = segments[2]

            if readprm == 'bond':
                type1 = segments[0]
                type2 = segments[1]
                b0 = float(segments[3])
                Kb = float(segments[2])
                b0 = float(b0)/10        # conversion from A -> nm
                Kb = Kb*2*kcal2kJ*1000/10    # converstion from kcal/mole/A**2 -> kJ/mole/nm**2 incl factor 2 (see definitions)
                type2bond[type1, type2] = [b0, Kb]

            if readprm == 'angle':
                type1 = segments[0]
                type2 = segments[1]
                type3 = segments[2]
                th0 = float(segments[4])
                cth = float(segments[3])
                cth = cth*2*kcal2kJ        # kJ/mol and an factor 2 (see definitions)
                if len(segments)>6:        # check for Urey-Bradley parameters
                    try:
                        S0  = float(segments[6])
                        Kub = float(segments[5])
                        S0  = S0/10
                        Kub = Kub*2*kcal2kJ*1000/10
                        ubFlag = True
                    except ValueError:
                        ubFlag = False
                else:
                    ubFlag = False
                if not ubFlag:
                    S0  = 0.0
                    Kub = 0.0
                type2angle[type1, type2, type3] = [th0, cth, S0, Kub]

            if readprm == 'dihedral':
                type1 = segments[0]
                type2 = segments[1]
                type3 = segments[2]
                type4 = segments[3]
                phi0 = float(segments[6])
                cp   = float(segments[4])
                mult = int(segments[5])
                cp   = cp*kcal2kJ            # conversion to kJ
                if type1 == 'X' and type4 == 'X':    # look for wildcards in positions 1 and 4
                    dihWilds.append([type2, type3, phi0, cp, mult])    # save them in a list
                else:                    # no wildcard - write to bon file
                    try:
                        type2dihedral[type1, type2, type3, type4] += [phi0, cp, mult]
                    except:
                        type2dihedral[type1, type2, type3, type4] = [phi0, cp, mult]
                    type2mult[type1, type2, type3, type4] = len(type2dihedral[type1, type2, type3, type4])/3

            if readprm == 'improper':
                type1 = segments[0]
                type2 = segments[1]
                type3 = segments[2]
                type4 = segments[3]
                q0 = float(segments [6])
                cq = float(segments [4])
                cq = cq*2*kcal2kJ        # conversion to kJ, factor 2 from definition difference
                if type2 == 'X' and type3 == 'X':    # look for wildcards in positions 2 and 3
                    impWilds.append([type1, type4, q0, cq])        # save them in a list
                else:                    # no wildcard - write to bon file
                    type2improper[type1, type2, type3, type4] = [q0, cq]

            if readprm == 'cmap':
                try:
                    segments[0] = float(segments[0])
                    for icmap in segments:
                        icmap = float(icmap)*kcal2kJ
                        cmapData[len(cmap)-1].append(icmap)
                except:
                    type1 = segments[0]
                    type2 = segments[1]
                    type3 = segments[2]
                    type4 = segments[3]
                    type5 = segments[7]
                    cmap.append([type1, type2, type3, type4, type5])
                    cmapCount.append(int(segments[8]))
                    cmapData.append([])


            if readprm == 'nonbonded':
                type     = segments[0]
                epsilon  = float(segments[2])
                RminHalf = float(segments[3])
                eps      = abs(epsilon*kcal2kJ)        # conversion to kJ and positive
                sigma    = 2*RminHalf/(10.0*2.0**(1.0/6.0))    # -> nm, double distance and rmin2sigma factor
                type2nonbonded[type] = [sigma, eps]
                if len(segments)> 6:        # test length to avoid IndexError
                    try:            # if possible, convert element 5 to float 
                        segments[5] = float(segments[5])
                    except:
                        None
                    # is segment 5 and 6 floats => there's 1-4 defined
                    if not isinstance(segments[5],str): # not string?
                        epsilon14  = float(segments[5])            # read charmm epsilon
                        eps14      = abs(epsilon14*kcal2kJ)        # conversion to gromacs units
                        Rmin14Half = float(segments[6])            # read charmm Rmin*1/2
                        sigma14    = 2*Rmin14Half/(10.0*2.0**(1.0/6.0))    # conversion to gromacs units
                        type2pair[type] = [sigma14, eps14]        # add to list

            if readprm == 'nbfix':
                type1 = segments[0]
                type2 = segments[1]
                epsilon = float(segments[2])
                Rmin    = float(segments[3])
                eps     = abs(epsilon*kcal2kJ)        # conversion to kJ and positive
                sigma   = Rmin/(10.0*2.0**(1.0/6.0))    # -> nm, double distance and rmin2sigma factor
                type2nbfix[type1,type2] = [sigma, eps]





#-------------------------------------------#
#               Read psfFile                #
#-------------------------------------------#

psfFile = open(sys.argv[2], 'r')

mol   = []
memb  = []
types = []

atoms     = []
bonds     = []
angles    = []
dihedrals = []
impropers = []
grps      = []
cmaps     = []

molinfo   = []
molNumber = {}

type2charge = {}

readpsf = None
for line in psfFile:
    line = line.upper().strip()
    if len(line) > 0:

        segments = line.split()

        if line.startswith('PSF'):
            xplor = False
            if 'XPLOR' in segments:
                xplor = True
            continue
        if line.endswith('NATOM'):
            readpsf = 'atom'
            continue
        if line.endswith('BONDS'):
            readpsf = 'bond'
            continue
        if line.endswith('ANGLES'):
            readpsf = 'angle'
            continue
        if line.endswith('DIHEDRALS'):
            readpsf = 'dihedral'
            continue
        if line.endswith('IMPROPERS'):
            readpsf = 'improper'
            continue
        if line.endswith('DONORS'):
            readpsf = None
            continue
        if line.endswith('ACCEPTORS'):
            readpsf = None
            continue
        if line.endswith('NNB'):
            readpsf = None
            continue
        if line.endswith('NST2'):
            readpsf = 'group'
            continue
        if line.endswith('MOLNT'):
            readpsf = None
            continue
        if line.endswith('NUMLPH'):
            readpsf = None
            continue
        if line.endswith('CROSS-TERMS'):
            readpsf = 'cmap'
            continue
        if readpsf == None:
            continue

        if readpsf == 'atom':
            atomid   = int(segments[0])
            segid    = segments[1]
            resid    = segments[2]
            resname  = segments[3]
            atomname = segments[4]
            atomtype = segments[5]
            atomchrg = float(segments[6])
            atommass = float(segments[7])
            if xplor:
                type = atomtype
            else:
                type = index2type[atomtype]
            if segid == resname or resname == 'TIP3' or segid == 'MEMB':
                try:
                    imol = mol.index(resname)
                except:
                    mol.append(resname)
                    molinfo.append({'firstres': resid, 'firstatom': atomname, 'firstatomid': atomid, 'natom': 0})
                    molNumber[resname] = 1

                    if segid == 'MEMB':
                        memb.append(resname)

                    atoms.append([])
                    bonds.append([])
                    angles.append([])
                    dihedrals.append([])
                    impropers.append([])
                    grps.append([])
                    cmaps.append([])

                    imol = len(mol) - 1

                if ( resid != molinfo[imol]['firstres'] or molNumber[resname] != 1 ) and atomname == molinfo[imol]['firstatom']:
                    molNumber[resname] += 1

                if molNumber[resname] == 1: 
                    atoms[imol].append({'type': type, 'resnr': resid, 'residu': resname, 'atom': atomname, 'charge': atomchrg, 'mass': atommass})
                    bonds[imol].append([])
                    molinfo[imol]['natom'] += 1

            else:
                try:
                    imol = mol.index(segid)
                except:
                    mol.append(segid)
                    molinfo.append({'firstatomid': atomid, 'natom': 0})
                    molNumber[segid] = 1

                    atoms.append([])
                    bonds.append([])
                    angles.append([])
                    dihedrals.append([])
                    impropers.append([])
                    grps.append([])
                    cmaps.append([])

                    imol = len(mol) - 1
                    
                atoms[imol].append({'type': type, 'resnr': resid, 'residu': resname, 'atom': atomname, 'charge': atomchrg, 'mass': atommass})
                bonds[imol].append([])
                molinfo[imol]['natom'] += 1

            try:
                types.index(type)
            except:
                types.append(type)
                type2charge[type] = atomchrg

        if readpsf == 'bond':
            bondNumber = len(segments)/2
            for i in range(bondNumber):
                atom1 = int(segments[2*i+0])
                atom2 = int(segments[2*i+1])
                for imol in range(len(mol)):
                    fatom = molinfo[imol]['firstatomid']
                    natom = molinfo[imol]['natom']
                    latom = fatom + natom - 1
                    if fatom <= atom1 <= latom and fatom <= atom2 <= latom:
                        atom1  -= fatom
                        atom2  -= fatom
                        bonds[imol][atom1].append(atom2)
                        bonds[imol][atom2].append(atom1)
                        break

        if readpsf == 'angle':
            angleNumber = len(segments)/3
            for i in range(angleNumber):
                atom1 = int(segments[3*i+0])
                atom2 = int(segments[3*i+1])
                atom3 = int(segments[3*i+2])
                for imol in range(len(mol)):
                    fatom = molinfo[imol]['firstatomid']
                    natom = molinfo[imol]['natom']
                    latom = fatom + natom - 1
                    if fatom <= atom1 <= latom and fatom <= atom2 <= latom and fatom <= atom3 <= latom:
                        atom1  -= fatom - 1
                        atom2  -= fatom - 1
                        atom3  -= fatom - 1
                        angles[imol].append([atom2, atom1, atom3])
                        break

        if readpsf == 'dihedral':
            dihedralNumber = len(segments)/4
            for i in range(dihedralNumber):
                atom1 = int(segments[4*i+0])
                atom2 = int(segments[4*i+1])
                atom3 = int(segments[4*i+2])
                atom4 = int(segments[4*i+3])
                for imol in range(len(mol)):
                    fatom = molinfo[imol]['firstatomid']
                    natom = molinfo[imol]['natom']
                    latom = fatom + natom - 1
                    if fatom <= atom1 <= latom and fatom <= atom2 <= latom and fatom <= atom3 <= latom and fatom <= atom4 <= latom:
                        atom1  -= fatom - 1
                        atom2  -= fatom - 1
                        atom3  -= fatom - 1
                        atom4  -= fatom - 1
                        dihedrals[imol].append([atom2, atom3, atom1, atom4])
                        break

        if readpsf == 'improper':
            improperNumber = len(segments)/4
            for i in range(improperNumber):
                atom1 = int(segments[4*i+0])
                atom2 = int(segments[4*i+1])
                atom3 = int(segments[4*i+2])
                atom4 = int(segments[4*i+3])
                for imol in range(len(mol)):
                    fatom = molinfo[imol]['firstatomid']
                    natom = molinfo[imol]['natom']
                    latom = fatom + natom - 1
                    if fatom <= atom1 <= latom and fatom <= atom2 <= latom and fatom <= atom3 <= latom and fatom <= atom4 <= latom:
                        atom1  -= fatom - 1
                        atom2  -= fatom - 1
                        atom3  -= fatom - 1
                        atom4  -= fatom - 1
                        impropers[imol].append([atom1, atom4, atom2, atom3])
                        break

        if readpsf == 'group':
            groupNumber = len(segments)/3
            for i in range(groupNumber):
                atom1 = int(segments[3*i+0]) + 1
                for imol in range(len(mol)):
                    fatom = molinfo[imol]['firstatomid']
                    natom = molinfo[imol]['natom']
                    latom = fatom + natom - 1
                    if fatom <= atom1 <= latom:
                        atom1  -= fatom - 1
                        grps[imol].append(atom1)
                        break

        if readpsf == 'cmap':
            atom1 = int(segments[0])
            atom2 = int(segments[1])
            atom3 = int(segments[2])
            atom4 = int(segments[3])
            atom5 = int(segments[7])
            for imol in range(len(mol)):
                fatom = molinfo[imol]['firstatomid']
                natom = molinfo[imol]['natom']
                latom = fatom + natom - 1
                if fatom <= atom1 <= latom and fatom <= atom2 <= latom and fatom <= atom3 <= latom and fatom <= atom4 <= latom and fatom <= atom5 <= latom:
                    atom1  -= fatom - 1
                    atom2  -= fatom - 1
                    atom3  -= fatom - 1
                    atom4  -= fatom - 1
                    atom5  -= fatom - 1
                    cmaps[imol].append([atom1, atom2, atom3, atom4, atom5])
                    break





#-------------------------------------------#
#               Write itpFiles              #
#-------------------------------------------#

# set the func parameter for bonds, angles and proper/improper dihedrals
# for further information see section 5.7.1 in gromacs documentation:
# ftp://ftp.gromacs.org/pub/manual/manual-4.6.5.pdf

funcForBonds     = 1
funcForAngles    = 5    # Urey-Bradley angle type
funcForDihedrals = 9    # special type for treating multiple entries (modification in source code)
funcForImpropers = 2
funcFor14        = 1    # 1-4 interaction pair type
funcForCmap      = 1

# particle type
ptype = 'A'

# dictionary for atom numbers, used for the nb file
element2atomNumber= {}
element2atomNumber['H']  = 1
element2atomNumber['HE'] = 2
element2atomNumber['LI'] = 3
element2atomNumber['B']  = 5
element2atomNumber['C']  = 6
element2atomNumber['N']  = 7
element2atomNumber['O']  = 8
element2atomNumber['F']  = 9
element2atomNumber['NE'] = 10
element2atomNumber['NA'] = 11
element2atomNumber['MG'] = 12
element2atomNumber['AL'] = 13
element2atomNumber['P']  = 15
element2atomNumber['S']  = 16
element2atomNumber['CL'] = 17
element2atomNumber['K']  = 19
element2atomNumber['CA'] = 20
element2atomNumber['FE'] = 26
element2atomNumber['CU'] = 29
element2atomNumber['ZN'] = 30
element2atomNumber['BR'] = 35
element2atomNumber['RB'] = 37
element2atomNumber['CD'] = 48
element2atomNumber['I']  = 53
element2atomNumber['CS'] = 55
element2atomNumber['BA'] = 56

mass2element = {}
mass2element['1.0'] = 'H'
mass2element['4.0'] = 'HE'
mass2element['6.9'] = 'LI'
mass2element['10.8'] = 'B'
mass2element['12.0'] = 'C'
mass2element['14.0'] = 'N'
mass2element['16.0'] = 'O'
mass2element['19.0'] = 'F'
mass2element['20.2'] = 'NE'
mass2element['23.0'] = 'NA'
mass2element['24.3'] = 'MG'
mass2element['27.0'] = 'AL'
mass2element['31.0'] = 'P'
mass2element['32.1'] = 'S'
mass2element['35.5'] = 'CL'
mass2element['39.1'] = 'K'
mass2element['40.1'] = 'CA'
mass2element['55.8'] = 'FE'
mass2element['63.5'] = 'CU'
mass2element['65.4'] = 'ZN'
mass2element['79.9'] = 'BR'
mass2element['85.5'] = 'RB'
mass2element['112.4'] = 'CD'
mass2element['126.9'] = 'I'
mass2element['132.9'] = 'CS'
mass2element['137.3'] = 'BA'

type2elementNumber = {}


# Build DataBase
#------------------------------

types.sort()

dbbonds     = []
dbpairs     = []
dbnbfix     = []
dbangles    = []
dbdihedrals = []
dbdihewilds = []
dbimpropers = []
dbimprwilds = []
dbcmap      = []

# bondtypes
for bond in type2bond:
    type1 = bond[0]
    type2 = bond[1]
    try:
        types.index(type1)
        types.index(type2)
        writeitp = True
    except:
        writeitp = False
    if writeitp:
        b0 = type2bond[type1,type2][0]
        Kb = type2bond[type1,type2][1]
        dbbonds.append([type1,type2,b0,Kb])

# pairtypes
for i in range(len(types)):
    for j in range(i,len(types)):
        type1 = types[i]
        type2 = types[j]
        ipair = False
        jpair = False
        nbfix = False
        try:
            type2nbfix[type1,type2]
            nbfix = True
        except:
            try:
                type2nbfix[type2,type1]
                nbfix = True
            except:
                nbfix = False
        if not nbfix:
            try:
                isigma14 = type2pair[type1][0]
                ieps14   = type2pair[type1][1]
                ipair    = True
            except:
                isigma14 = type2nonbonded[type1][0]
                ieps14   = type2nonbonded[type1][1]
            try:
                jsigma14 = type2pair[type2][0]
                jeps14   = type2pair[type2][1]
                jpair    = True
            except:
                jsigma14 = type2nonbonded[type2][0]
                jeps14   = type2nonbonded[type2][1]
            if ipair or jpair:
                sigma14 = (isigma14 + jsigma14)/2.0
                eps14   = (ieps14 * jeps14)**0.5
                dbpairs.append([type1,type2,sigma14,eps14])

# nbfix
for nbfix in type2nbfix:
    type1 = nbfix[0]
    type2 = nbfix[1]
    sigma = type2nbfix[type1,type2][0]
    eps   = type2nbfix[type1,type2][1]
    try:
        types.index(type1)
        types.index(type2)
        try:
            dbnbfix.index([type1,type2,sigma,eps])
        except:
            dbnbfix.append([type1,type2,sigma,eps])
    except:
        None

# angletypes
for angle in type2angle:
    type1 = angle[0]
    type2 = angle[1]
    type3 = angle[2]
    try:
        types.index(type1)
        types.index(type2)
        types.index(type3)
        writeitp = True
    except:
        writeitp = False
    if writeitp:
        th0 = type2angle[type1,type2,type3][0]
        cth = type2angle[type1,type2,type3][1]
        S0  = type2angle[type1,type2,type3][2]
        Kub = type2angle[type1,type2,type3][3]
        dbangles.append([type1,type2,type3,th0,cth,S0,Kub])

# dihedraltypes
for dihedral in type2dihedral:
    type1 = dihedral[0]
    type2 = dihedral[1]
    type3 = dihedral[2]
    type4 = dihedral[3]
    try:
        types.index(type1)
        types.index(type2)
        types.index(type3)
        types.index(type4)
        writeitp = True
    except:
        writeitp = False
    if writeitp:
        multNumber = type2mult[type1,type2,type3,type4]
        if multNumber >= 1:
            for m in range(multNumber):
                phi  = type2dihedral[type1,type2,type3,type4][m*3]
                cp   = type2dihedral[type1,type2,type3,type4][m*3+1]
                mult = type2dihedral[type1,type2,type3,type4][m*3+2]
                dbdihedrals.append([type1,type2,type3,type4,phi,cp,mult])
for wild in dihWilds:
    type2 = wild[0]
    type3 = wild[1]
    try:
        types.index(type2)
        types.index(type3)
        writeitp = True
    except:
        writeitp = False
    if writeitp:
        phi  = wild[2]
        cp   = wild[3]
        mult = wild[4]
        dbdihewilds.append(['X',type2,type3,'X',phi,cp,mult])

# impropertypes
for improper in type2improper:
    type1 = improper[0]
    type2 = improper[1]
    type3 = improper[2]
    type4 = improper[3]
    try:
        types.index(type1)
        types.index(type2)
        types.index(type3)
        types.index(type4)
        writeitp = True
    except:
        writeitp = False
    if writeitp:
        q0 = type2improper[type1,type2,type3,type4][0]
        cq = type2improper[type1,type2,type3,type4][1]
        dbimpropers.append([type1,type2,type3,type4,q0,cq])
for wild in impWilds:
    type1 = wild[0]
    type4 = wild[1]
    try:
        types.index(type1)
        types.index(type4)
        writeitp = True
    except:
        writeitp = False
    if writeitp:
        q0 = wild[2]
        cq = wild[3]
        dbimprwilds.append([type1,'X','X',type4,q0,cq])

# cmaptypes
for icmap in range(len(cmap)):
    type1 = cmap[icmap][0]
    type2 = cmap[icmap][1]
    type3 = cmap[icmap][2]
    type4 = cmap[icmap][3]
    type5 = cmap[icmap][4]
    try:
        types.index(type1)
        types.index(type2)
        types.index(type3)
        types.index(type4)
        types.index(type5)
        dbcmap.append(icmap)
    except:
        None

dbbonds.sort()
dbpairs.sort()
dbnbfix.sort()
dbangles.sort()
dbdihedrals.sort()
dbdihewilds.sort()
dbimpropers.sort()
dbimprwilds.sort()


# Write parFiles: charmm36.itp
#------------------------------
itpFile = open('gromacs/charmm36.itp', 'w')

itpFile.write(';;\n')
itpFile.write(';; Generated by CHARMM-GUI (http://www.charmm-gui.org) v1.7\n')
itpFile.write(';;\n')
itpFile.write(';; psf2itp_mol.py\n')
itpFile.write(';;\n')
itpFile.write(';; Correspondance:\n')
itpFile.write(';; jul316@lehigh.edu or wonpil@lehigh.edu\n')
itpFile.write(';;\n')
itpFile.write(';; CHARMM36 FF in GROMACS format\n')
itpFile.write(';;\n\n')

# defaults
itpFile.write('\n[ defaults ]\n')
itpFile.write('; nbfunc\tcomb-rule\tgen-pairs\tfudgeLJ\tfudgeQQ\n')
itpFile.write('1\t2\tyes\t1.0\t1.0\n')

# atomtypes
itpFile.write('\n[ atomtypes ]\n')
itpFile.write('; name\tat.num\tmass\tcharge\tptype\tsigma\tepsilon\t;\tsigma_14\tepsilon_14\n')
for type in types:
    mass = type2mass[type]
    charge = type2charge[type]
    element = mass2element[str(round(mass,1))]
    elementNumber = element2atomNumber[element]
    type2elementNumber[type] = elementNumber
    sigma = type2nonbonded[type][0]
    eps = type2nonbonded[type][1]
    try:
        sigma14 = type2pair[type][0]
        eps14 = type2pair[type][1]
        itpFile.write(' %7s %5d %10.4f %10.3f %5s %20.11e %15.6e ; %19.11e %15.6e \n' % (type, elementNumber, mass, charge, ptype, sigma, eps, sigma14, eps14))
    except:
        itpFile.write(' %7s %5d %10.4f %10.3f %5s %20.11e %15.6e \n' % (type, elementNumber, mass, charge, ptype, sigma, eps))

# nbfix
if len(dbnbfix) > 0:
    itpFile.write('\n[ nonbond_params ]\n')
    itpFile.write('; i\tj\tfunc\tsigma\tepsilon\n')
    for nbfix in dbnbfix:
        type1   = nbfix[0]
        type2   = nbfix[1]
        sigma   = nbfix[2]
        eps     = nbfix[3]
        itpFile.write('%7s %7s %5d %18.11e %18.11e \n' % (type1, type2, funcFor14, sigma, eps))

# bondtypes
if len(dbbonds) > 0:
    itpFile.write('\n[ bondtypes ]\n')
    itpFile.write('; i\tj\tfunc\tb0\tKb\n')
    for bond in dbbonds:
        type1 = bond[0]
        type2 = bond[1]
        b0    = bond[2] 
        Kb    = bond[3]
        itpFile.write('%7s %7s %5d %13.6e %13.6e\n' % (type1, type2, funcForBonds, b0, Kb))

# pairtypes
if len(dbpairs) > 0:
    itpFile.write('\n[ pairtypes ]\n')
    itpFile.write('; i\tj\tfunc\tsigma1-4\tepsilon1-4\n')
    for pair in dbpairs:
        type1   = pair[0]
        type2   = pair[1]
        sigma14 = pair[2]
        eps14   = pair[3]
        itpFile.write('%7s %7s %5d %18.11e %18.11e \n' % (type1, type2, funcFor14, sigma14, eps14))

# angletypes
if len(dbangles) > 0:
    itpFile.write('\n[ angletypes ]\n')
    itpFile.write('; i\tj\tk\tfunc\tth0\tcth\tS0\tKub\n')
    for angle in dbangles:
        type1 = angle[0]
        type2 = angle[1]
        type3 = angle[2]
        th0   = angle[3]
        cth   = angle[4]
        S0    = angle[5]
        Kub   = angle[6] 
        itpFile.write('%7s %7s %7s %5d %14.7e %14.7e %14.7e %14.7e\n' % (type1, type2, type3, funcForAngles, th0, cth, S0, Kub))

# dihedraltypes
if len(dbdihedrals) > 0 or len(dbdihewilds) > 0:
    itpFile.write('\n[ dihedraltypes ]\n')
    itpFile.write('; i\tj\tk\tl\tfunc\tphi0\tcp\tmult\n')
    for dihedral in dbdihedrals:
        type1 = dihedral[0]
        type2 = dihedral[1]
        type3 = dihedral[2]
        type4 = dihedral[3]
        phi   = dihedral[4]
        cp    = dihedral[5]
        mult  = dihedral[6]
        itpFile.write('%7s %7s %7s %7s %5d %13.6e %13.6e %6d\n' % (type1, type2, type3, type4, funcForDihedrals, phi, cp, mult))
    for dihedral in dbdihewilds:
        type1 = dihedral[0]
        type2 = dihedral[1]
        type3 = dihedral[2]
        type4 = dihedral[3]
        phi   = dihedral[4]
        cp    = dihedral[5]
        mult  = dihedral[6]
        itpFile.write('%7s %7s %7s %7s %5d %13.6e %13.6e %6d\n' % (type1, type2, type3, type4, funcForDihedrals, phi, cp, mult))

# impropertypes
if len(dbimpropers) > 0 or len(dbimprwilds) > 0:
    itpFile.write('\n[ dihedraltypes ]\n')
    itpFile.write('; i\tj\tk\tl\tfunc\tq0\tcq\n')
    for improper in dbimpropers:
        type1 = improper[0]
        type2 = improper[1]
        type3 = improper[2]
        type4 = improper[3]
        q0    = improper[4]
        cq    = improper[5]
        itpFile.write('%7s %7s %7s %7s %5d %13.6e %13.6e\n' % (type1, type2, type3, type4, funcForImpropers, q0, cq))
    for improper in dbimprwilds:
        type1 = improper[0]
        type2 = improper[1]
        type3 = improper[2]
        type4 = improper[3]
        q0    = improper[4]
        cq    = improper[5]
        itpFile.write('%7s %7s %7s %7s %5d %13.6e %13.6e\n' % (type1, type2, type3, type4, funcForImpropers, q0, cq))

# cmaptypes
if len(dbcmap) > 0:
    itpFile.write('\n[ cmaptypes ]\n')
    itpFile.write('; i j k l m\n')
    for icmap in dbcmap:
        type1 = cmap[icmap][0]
        type2 = cmap[icmap][1]
        type3 = cmap[icmap][2]
        type4 = cmap[icmap][3]
        type5 = cmap[icmap][4]
        ncmap = cmapCount[icmap]
        itpFile.write('%s %s %s %s %s %d %d %d\\\n' % (type1, type2, type3, type4, type5, funcForCmap, ncmap, ncmap))
        for icmapdata in range(len(cmapData[icmap])):
            itpFile.write('%f' % cmapData[icmap][icmapdata])
            if (icmapdata+1) == len(cmapData[icmap]):
                itpFile.write('\n\n')
            else:
                if (icmapdata+1) % 10 == 0:
                    itpFile.write('\\\n')
                else:
                    itpFile.write(' ')

itpFile.close()


# Write itpFiles
#--------------------

# molecule.itp
for imol, molname in enumerate(mol):
    molatoms     = atoms[imol]
    molbonds     = bonds[imol]
    molangles    = angles[imol]
    moldihedrals = dihedrals[imol]
    molimpropers = impropers[imol]
    molgrps      = grps[imol]
    molcmaps     = cmaps[imol]

    # build 1-4 pairs
    molpairs = []
    for i in range(len(molbonds)):
        for j in molbonds[i]:
            for k in molbonds[j]:
                if k != i:
                    for l in molbonds[k]:
                        if l != j and l > i:
                            pairpass = False
                            for j2 in molbonds[i]:
                                if l in molbonds[j2]:
                                    pairpass = True
                            if not pairpass:
                                try:
                                    molpairs.index([i, l])
                                except:
                                    molpairs.append([i, l])

    itpFile = open('gromacs/'+molname+'.itp', 'w')

    nrexcl = 1
    if len(molbonds) > 1:
        nrexcl = 2
    if len(moldihedrals) > 0:
        nrexcl = 3

    itpFile.write(';;\n')
    itpFile.write(';; Generated by CHARMM-GUI (http://www.charmm-gui.org) v1.7\n')
    itpFile.write(';;\n')
    itpFile.write(';; psf2itp_mol.py\n')
    itpFile.write(';;\n')
    itpFile.write(';; Correspondance:\n')
    itpFile.write(';; jul316@lehigh.edu or wonpil@lehigh.edu\n')
    itpFile.write(';;\n')
    itpFile.write(';; GROMACS topology file for %s\n' % molname.upper())
    itpFile.write(';;\n\n')

    itpFile.write('\n[ moleculetype ]\n')
    itpFile.write('; name\tnrexcl\n')
    itpFile.write('%s\t %5d\n' % (molname.upper(), nrexcl))

    # atoms
    itpFile.write('\n[ atoms ]\n')
    itpFile.write('; nr\ttype\tresnr\tresidu\tatom\tcgnr\tcharge\tmass\n')
    qtot = 0.0
    for i,atom in enumerate(molatoms):
        type   = atom['type']
        resnr  = atom['resnr']
        residu = atom['residu']
        iatom  = atom['atom']
        charge = atom['charge']
        qtot  += charge
        try:
            group = molgrps.index(i+1) + 1
        except:
            None
        itpFile.write(' %5d %10s %6s %8s %6s %6d %10.3f %10.4f   ; qtot %6.3f\n' % (i+1, type, resnr, residu, iatom, i+1, charge, type2mass[type], qtot))

    if molname != 'TIP3':

        # bonds
        if len(molbonds[0]) > 0:
            itpFile.write('\n[ bonds ]\n')
            itpFile.write('; ai\taj\tfunct\tb0\tKb\n')
            for i, bond in enumerate(molbonds):
                bond.sort()
                for j in bond:
                    if j > i:
                        atom1 = i+1
                        atom2 = j+1
                        itpFile.write('%5d %5d %5d\n' % (atom1, atom2, funcForBonds))

        # pairs
        if len(molpairs) > 0:
            itpFile.write('\n[ pairs ]\n')
            itpFile.write('; ai\taj\tfunct\tc6\tc12\n')
            molpairs.sort()
            for i in range(len(molpairs)):
                atom1 = molpairs[i][0]
                atom4 = molpairs[i][1]
                itpFile.write('%5d %5d %5d \n' % (atom1+1, atom4+1, funcFor14))

        # angles
        if len(molangles) > 0:
            itpFile.write('\n[ angles ] \n')
            itpFile.write('; ai\taj\tak\tfunct\tth0\tcth\tS0\tKub\n')
            molangles.sort()
            for angle in molangles:
                atom1 = angle[1]
                atom2 = angle[0]
                atom3 = angle[2]
                itpFile.write('%5d %5d %5d %5d\n' % (atom1, atom2, atom3, funcForAngles))

        # dihedrals
        if len(moldihedrals) > 0:
            itpFile.write('\n[ dihedrals ]\n')
            itpFile.write('; ai\taj\tak\tal\tfunct\tphi0\tcp\tmult\n')
            moldihedrals.sort()
            for dihedral in moldihedrals:
                atom1 = dihedral[2]
                atom2 = dihedral[0]
                atom3 = dihedral[1]
                atom4 = dihedral[3]
                itpFile.write('%5d %5d %5d %5d %5d\n' % (atom1, atom2, atom3, atom4, funcForDihedrals))

        # impropers
        if len(molimpropers) > 0:
            itpFile.write('\n[ dihedrals ]\n')
            itpFile.write('; ai\taj\tak\tal\tfunct\tq0\tcq\n')
            molimpropers.sort()
            for improper in molimpropers:
                atom1 = improper[0]
                atom2 = improper[2]
                atom3 = improper[3]
                atom4 = improper[1]
                itpFile.write('%5d %5d %5d %5d %5d\n' % (atom1, atom2, atom3, atom4, funcForImpropers))

        # cmaps
        if len(molcmaps) > 0:
            itpFile.write('\n[ cmap ] \n')
            itpFile.write('; ai\taj\tak\tal\tam\tfunct\n')
            molcmaps.sort()
            for icmap in molcmaps:
                atom1 = icmap[0]
                atom2 = icmap[1]
                atom3 = icmap[2]
                atom4 = icmap[3]
                atom5 = icmap[4]
                itpFile.write('%5d %5d %5d %5d %5d %5d\n' % (atom1, atom2, atom3, atom4, atom5, funcForCmap))

        if molname != residu or molname in memb:
           # position or dihedral restraints for proteins, carbohydrates, and ligands
           itpFile.write('\n#ifdef REST_ON\n')
           itpFile.write('#include "%s_rest.itp"\n' % molname.upper())
           itpFile.write('#endif\n')

    else:
        try:
            b0_1 = type2bond['OT','HT'][0]
            b0_2 = type2bond['HT','HT'][0]
        except:
            b0_1 = type2bond['HT','OT'][0]
            b0_2 = type2bond['HT','HT'][0]

        itpFile.write('\n[ settles ]\n')
        itpFile.write('; i\tj\tfunct\tlength\n')
        itpFile.write('1\t1\t%13.6e\t%13.6e\n' % (b0_1, b0_2))

        itpFile.write('\n[ exclusions ]\n')
        itpFile.write('1 2 3\n')
        itpFile.write('2 1 3\n')
        itpFile.write('3 1 2\n')

    itpFile.close()


# Write topFile
#--------------------

# topol.top

topFile = open('gromacs/topol.top', 'w')

topFile.write(';;\n')
topFile.write(';; Generated by CHARMM-GUI (http://www.charmm-gui.org) v1.7\n')
topFile.write(';;\n')
topFile.write(';; psf2itp_mol.py\n')
topFile.write(';;\n')
topFile.write(';; Correspondance:\n')
topFile.write(';; jul316@lehigh.edu or wonpil@lehigh.edu\n')
topFile.write(';;\n')
topFile.write(';; The main GROMACS topology file\n')
topFile.write(';;\n\n')

topFile.write('; Include forcefield parameters\n')
topFile.write('#include "charmm36.itp"\n')
for imol in mol:
    topFile.write('#include "%s.itp"\n' % imol.upper())

topFile.write('\n[ system ]\n')
topFile.write('; Name\n')
topFile.write('Title\n')

topFile.write('\n[ molecules ]\n')
topFile.write('; Compound\t#mols\n')
for imol in mol:
    topFile.write('%-6s\t%12d\n' % (imol,molNumber[imol]))

topFile.close()

