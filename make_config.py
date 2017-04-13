import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-nx", type=int, default=20, help="number of ligands in x-direction")
parser.add_argument("-nz", type=int, default=12, help="number of ligands in z-direction")
parser.add_argument("-fmt", type=str, choices=["xyz", "lmp"], default='xyz', help="output format")
args = parser.parse_args()

data = np.genfromtxt('/home/pratima/Nanoparticle-Rotation/diatomic-probe/backbone.data', delimiter=' ')
bottom = np.genfromtxt('/home/pratima/Nanoparticle-Rotation/diatomic-probe/lowerlig.data', delimiter=' ')
# liq_hex = np.genfromtxt('/home/pratima/OleicAcidLigand/liq_solv.txt', delimiter=' ')
liq_hex = np.genfromtxt('/home/pratima/Nanoparticle-Rotation/diatomic-probe/hex_solv.data', delimiter=' ')
dye = np.genfromtxt('/home/pratima/Nanoparticle-Rotation/diatomic-probe/dye.data', delimiter=' ')
x0 = data[0,1]
y0 = data[0,2]
z0 = data[0,3]
bx0 = bottom[0,1]
by0 = bottom[0,2]
bz0 = bottom[0,3]
sx0 = np.min(liq_hex[:,1])
sy0 = np.min(liq_hex[:,2])
sz0 = np.min(dye[:,3])
dx0 = np.min(dye[:,1])
dy0 = np.min(dye[:,2])
dz0 = np.min(dye[:,3])
centred_data = np.zeros( (len(data), 3) )
centred_bottom = np.zeros( (len(bottom), 3) )
centred_hex = np.zeros( (len(liq_hex), 3) )
centred_dye = np.zeros( (len(dye), 3) )
nx = args.nx
nz = args.nz
n_molecules = nx * nz
n_bottom = len(bottom)
n_chain = 18
n_bottommolecules = n_bottom / n_chain
n_bonds = n_chain - 1
n_angles = n_chain - 2
n_dihedrals = n_chain - 3
n_hexatoms = len(liq_hex) 
n_hex = 6
n_hexmolecules = n_hexatoms / n_hex
n_hexbonds = n_hex - 1
n_hexangles = n_hex - 2
n_hexdihedrals = n_hex - 3
n_dye = len(dye)
n_dyebonds = n_dye + 2
n_dyeangles = n_dye + 6
n_dyedihedrals = 9
n_dyeimpropers = 5
xlo =  -5.0
xhi =  80.0
ylo = -124.0
yhi =  100.0
zlo =  -5.0
zhi =  80.0

# print n_hexdihedrals * n_hexmolecules
# print n_bottommolecules * n_dihedrals

rotate = np.array([ [np.cos(0.5 * np.pi), 0, -np.sin(0.5 * np.pi)], [0, 1, 0], [np.sin(0.5 * np.pi), 0,  np.cos(0.5 * np.pi)] ])

# centre the original backbone
for i in range(len(data)):
    centred_data[i,0] = data[i,1] - x0
    centred_data[i,1] = data[i,2] - y0
    centred_data[i,2] = data[i,3] - z0

# centre the original backbone
for i in range(len(bottom)):
    centred_bottom[i,0] = bottom[i,1] - bx0
    centred_bottom[i,1] = bottom[i,2] - by0
    centred_bottom[i,2] = bottom[i,3] - bz0

# centre the hex chain
for i in range(len(liq_hex)):
    centred_hex[i,0] = liq_hex[i,1]
    centred_hex[i,1] = liq_hex[i,2] + 20.0
    centred_hex[i,2] = liq_hex[i,3]

# centre the hex chain
for i in range(len(dye)):
    centred_dye[i,0] = dye[i,1] - dx0 + 4.12 * 10
    centred_dye[i,1] = dye[i,2] - dy0
    centred_dye[i,2] = dye[i,3] - dz0 + 6.75 * 6

# rotate the dye molecule to be upright
# for i in range(len(centred_dye)):
#     centred_dye[i, :] = rotate.dot(centred_dye[i, :])

# print sy0, min(centred_hex[:,1])

if args.fmt == "xyz":
    print "{}".format( (n_molecules - 5) * n_chain + n_molecules * 16 + 2 * n_hexatoms + n_bottom + n_dye)
    print "Comment Line"
    nmol = 0

    # dye on surface
    for i in range(n_dye):
        if (i == 6):
            print "3	{:4.5f}	{:4.5f}	{:4.5f}".format( centred_dye[i,0], centred_dye[i,1], centred_dye[1,2] )
        elif (i == 7): 
            print "4	{:4.5f}	{:4.5f}	{:4.5f}".format( centred_dye[i,0], centred_dye[i,1], centred_dye[1,2] )
        elif (i == 8):
            print "5	{:4.5f}	{:4.5f}	{:4.5f}".format( centred_dye[i,0], centred_dye[i,1], centred_dye[1,2] )
        elif (i == 12 or i == 13 or i ==9 or i == 10):
            print "6	{:4.5f}	{:4.5f}	{:4.5f}".format( centred_dye[i,0], centred_dye[i,1], centred_dye[1,2] )
        elif (i > 8 and i < 22):
            print "7	{:4.5f}	{:4.5f}	{:4.5f}".format( centred_dye[i,0], centred_dye[i,1], centred_dye[1,2] )
        else:
            print "2	{:4.5f}	{:4.5f}	{:4.5f}".format( centred_dye[i,0], centred_dye[i,1], centred_dye[1,2] )

    # atomic coordinates
    for i in range(nx):
        for j in range(nz):
            coords = np.copy(centred_data)
            # shift molecule
            coords[:, 0] = coords[:, 0] + 4.12 * i
            coords[:, 2] = coords[:, 2] + 6.75 * j

            # skip dye position
            if (i >= 10 and i <=14 and j == 6):
                continue
       
            for k in range(n_chain):
                if (k == n_chain-1):
                    print "1       {:4.5f} {:4.5f} {:4.5f}".format( coords[k,0], coords[k,1], coords[k,2] )
                else:
                    print "2       {:4.5f} {:4.5f} {:4.5f}".format( coords[k,0], coords[k,1], coords[k,2] )
    
            nmol = nmol + 1

    # Cd-S surface coordinates
    for i in range(nx):
        for j in range(nz):
            coords = centred_data[0, :] - np.array([0.0, 1.56, 0.0])
            coords[0] = coords[0] + 4.12 * i
            coords[2] = coords[2] + 6.75 * j
            print "9	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0], coords[1], coords[2] )
            print "8	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0], coords[1], coords[2] + 2.5299 )
            print "9	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0] - 2.06073, coords[1] - 1.18977, coords[2] + 3.375 )
            print "8	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0] - 2.06073, coords[1] - 1.18977, coords[2] - 3.375 + 2.5299 )
            print "9	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0] - 2.06073, coords[1] - 3.5693, coords[2] )
            print "8	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0] - 2.06073, coords[1] - 3.5693, coords[2] + 2.5299 )
            print "9	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0], coords[1] - 3.5693 - 1.18977, coords[2] + 3.375 )
            print "8	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0], coords[1] - 3.5693 - 1.18977, coords[2] - 3.375 + 2.5299 )

            print "9	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0], coords[1] - 7.13848, coords[2] )
            print "8	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0], coords[1] - 7.13848, coords[2] + 2.5299 )
            print "9	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0] - 2.06073, coords[1] - 7.13848 - 1.18977, coords[2] + 3.375 )
            print "8	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0] - 2.06073, coords[1] - 7.13848 - 1.18977, coords[2] - 3.375 + 2.5299 )
            print "9	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0] - 2.06073, coords[1] - 7.13848 - 3.5693, coords[2] )
            print "8	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0] - 2.06073, coords[1] - 7.13848 - 3.5693, coords[2] + 2.5299 )
            print "9	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0], coords[1] - 7.13848 - 3.5693 - 1.18977, coords[2] + 3.375 )
            print "8	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[0], coords[1] - 7.13848 - 3.5693 - 1.18977, coords[2] - 3.375 + 2.5299 )

    # hex molecules here
    coords = centred_hex + np.array([0.0, 0.0, 0.0])
    for i in range(n_hexatoms):    
        if ((i%6 == 0) or ((i+1)%6 == 0)):
            print "10	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[i,0], coords[i,1], coords[i,2] )
        else:
            print "11	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[i,0], coords[i,1], coords[i,2] )

    # bottom oleic acid ligands here
    coords = centred_bottom + np.array([0.0, 0.0, 0.0])
    for i in range(n_bottom):
        if ((i% (n_chain-1) ) == 0):
            print "1       {:4.5f} {:4.5f} {:4.5f}".format( coords[i,0], coords[i,1] + -15.457550, coords[i,2] )
        else:
            print "2       {:4.5f} {:4.5f} {:4.5f}".format( coords[i,0], coords[i,1] + -15.457550, coords[i,2] )

    # bottom hex molecules here
    coords = centred_hex + np.array([0.0, -55.5 - 80, 0.0])
    for i in range(n_hexatoms):    
        if ((i%6 == 0) or ((i+1)%6 == 0)):
            print "10	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[i,0], coords[i,1], coords[i,2] )
        else:
            print "11	{:4.5f}	{:4.5f}	{:4.5f}".format( coords[i,0], coords[i,1], coords[i,2] )

if args.fmt == "lmp":
    # beginning of LAMMPS input file
    print "LAMMPS coords for {} oleic acid-like molecules\n".format( n_molecules )
    print "{}	atoms".format( (n_molecules - 5) * n_chain + n_molecules * 16 + 2 * n_hexatoms + n_bottom + n_dye )
    print "{}	bonds".format( n_molecules * n_bonds + n_molecules + 2 * n_hexmolecules * n_hexbonds + n_bottommolecules * n_bonds+ n_dyebonds )
    print "{}	angles".format( n_molecules * n_angles + 2 * n_hexmolecules * n_hexangles + n_bottommolecules * n_angles + n_dyeangles )
    print "{}	dihedrals".format( n_molecules * n_dihedrals + 2 * n_hexmolecules * n_hexdihedrals + n_bottommolecules * n_dihedrals + n_dyedihedrals )
    print "{} impropers".format( n_dyeimpropers )
    
    print "11	atom types"
    print "3	bond types"
    print "3	angle types"
    print "3	dihedral types"
    print "1	improper types"
    
    print "{} {}	xlo xhi".format(xlo, xhi)
    print "{} {}	ylo yhi".format(ylo, yhi)
    print "{} {}	zlo zhi\n".format(zlo, zhi)
    
    print "Masses\n"
    print "1 15.034"
    print "2 14.026"
    print "3 112.411"
    print "4 32.065"
    print "5 15.034"
    print "6 14.026\n"
    print "Atoms\n"

    nmol = 0
    # begin defining actual coordinates
    for i in range(nx):
        for j in range(nz):
            coords = np.copy(centred_data)
            # shift molecule
            coords[:, 0] = coords[:, 0] + 4.12 * i
            coords[:, 2] = coords[:, 2] + 6.75 * j
        
            for k in range(n_chain):
                if (k == n_chain-1):
                    print "{}	{}	1	{:4.5f}	{:4.5f}	{:4.5f}".format( nmol * n_chain + k + 1, nmol+1, coords[k,0], coords[k,1], coords[k,2]  )
                else:
                    print "{}	{}	2	{:4.5f}	{:4.5f}	{:4.5f}".format( nmol * n_chain + k + 1, nmol+1, coords[k,0], coords[k,1], coords[k,2]  )
            nmol = nmol + 1

    # Cd-S surface coordinates
    ncd = 0
    for i in range(nx):
        for j in range(nz):
            coords = centred_data[0, :] - np.array([0.0, 1.56, 0.0])
            coords[0] = coords[0] + 4.12 * i
            coords[2] = coords[2] + 6.75 * j
            print "{}	{}	3	{:4.5f}	{:4.5f}	{:4.5f}".format( nmol * n_chain + ncd + 1, nmol+1, coords[0], coords[1], coords[2] )
            print "{}	{}	4	{:4.5f}	{:4.5f}	{:4.5f}".format( nmol * n_chain + ncd + 2, nmol+1, coords[0], coords[1], coords[2] + 2.5299 )
            ncd = ncd + 2

    for i in range(nx):
        for j in range(nz):
            coords = centred_data[0, :] - np.array([0.0, 1.56, 0.0])
            coords[0] = coords[0] + 4.12 * i
            coords[2] = coords[2] + 6.75 * j
            print "{}	{}	3	{:4.5f}	{:4.5f}	{:4.5f}".format( nmol * n_chain + ncd + 1, nmol+1, coords[0] - 2.06073, coords[1] - 1.18977, coords[2] + 3.375 )
            print "{}	{}	4	{:4.5f}	{:4.5f}	{:4.5f}".format( nmol * n_chain + ncd + 2, nmol+1, coords[0] - 2.06073, coords[1] - 1.18977, coords[2] - 3.375 + 2.5299 )
            print "{}	{}	3	{:4.5f}	{:4.5f}	{:4.5f}".format( nmol * n_chain + ncd + 3, nmol+1, coords[0] - 2.06073, coords[1] - 3.5693, coords[2] )
            print "{}	{}	4	{:4.5f}	{:4.5f}	{:4.5f}".format( nmol * n_chain + ncd + 4, nmol+1, coords[0] - 2.06073, coords[1] - 3.5693, coords[2] + 2.5299 )
            print "{}	{}	3	{:4.5f}	{:4.5f}	{:4.5f}".format( nmol * n_chain + ncd + 5, nmol+1, coords[0], coords[1] - 3.5693 - 1.18977, coords[2] + 3.375 )
            print "{}	{}	4	{:4.5f}	{:4.5f}	{:4.5f}".format( nmol * n_chain + ncd + 6, nmol+1, coords[0], coords[1] - 3.5693 - 1.18977, coords[2] - 3.375 + 2.5299 )

            print "{}	{}	3	{:4.5f}	{:4.5f}	{:4.5f}".format( nmol * n_chain + ncd + 7, nmol+1, coords[0], coords[1] - 7.13848, coords[2] )
            print "{}	{}	4	{:4.5f}	{:4.5f}	{:4.5f}".format( nmol * n_chain + ncd + 8, nmol+1, coords[0], coords[1] - 7.13848, coords[2] + 2.5299 )
            print "{}	{}	3	{:4.5f}	{:4.5f}	{:4.5f}".format( nmol * n_chain + ncd + 9, nmol+1, coords[0] - 2.06073, coords[1] - 7.13848 - 1.18977, coords[2] + 3.375 )
            print "{}	{}	4	{:4.5f}	{:4.5f}	{:4.5f}".format( nmol * n_chain + ncd + 10, nmol+1, coords[0] - 2.06073, coords[1] - 7.13848 - 1.18977, coords[2] - 3.375 + 2.5299 )
            print "{}	{}	3	{:4.5f}	{:4.5f}	{:4.5f}".format( nmol * n_chain + ncd + 11, nmol+1, coords[0] - 2.06073, coords[1] - 7.13848 - 3.5693, coords[2] )
            print "{}	{}	4	{:4.5f}	{:4.5f}	{:4.5f}".format( nmol * n_chain + ncd + 12, nmol+1, coords[0] - 2.06073, coords[1] - 7.13848 - 3.5693, coords[2] + 2.5299 )
            print "{}	{}	3	{:4.5f}	{:4.5f}	{:4.5f}".format( nmol * n_chain + ncd + 13, nmol+1, coords[0], coords[1] - 7.13848 - 3.5693 - 1.18977, coords[2] + 3.375 )
            print "{}	{}	4	{:4.5f}	{:4.5f}	{:4.5f}".format( nmol * n_chain + ncd + 14, nmol+1, coords[0], coords[1] - 7.13848 - 3.5693 - 1.18977, coords[2] - 3.375 + 2.5299 )
            ncd = ncd + 14

    # hex molecules here
    num_mol = nmol + 2
    num_atom = n_molecules * n_chain + n_molecules * 16 + 1

    coords = centred_hex + np.array([0.0, 0.0, 0.0]) 
    for i in range(n_hexatoms):
        if ((i%6 == 0) or ((i+1)%6 == 0)):
            print "{}	{}	5	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[i,0], coords[i,1], coords[i,2] )
        else:
            print "{}	{}	6	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[i,0], coords[i,1], coords[i,2] )

        num_atom = num_atom + 1
        # increase molecule number every 6 atoms
        if ((i+1)%6 == 0):
            num_mol = num_mol + 1

    # bottom molecules here
    coords = centred_bottom + np.array([0.0, 0.0, 0.0])
    for i in range(n_bottom):
        if ((i% (n_chain-1) ) == 0):
            print "{}	{}	1       {:4.5f} {:4.5f} {:4.5f}".format( num_atom, num_mol, coords[i,0], coords[i,1] + -15.457550, coords[i,2] )
        else:
            print "{}	{}	2       {:4.5f} {:4.5f} {:4.5f}".format( num_atom, num_mol, coords[i,0], coords[i,1] + -15.457550, coords[i,2] )
        num_atom = num_atom + 1
        if (i%n_chain == 0):
            num_mol = num_mol + 1

    # bottom hex molecules here
    coords = centred_hex + np.array([0.0, -55.5 - 80, 0.0])
    for i in range(n_hexatoms):    
        if ((i%6 == 0) or ((i+1)%6 == 0)):
            print "{}	{}	5	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[i,0], coords[i,1], coords[i,2] )
        else:
            print "{}	{}	6	{:4.5f}	{:4.5f}	{:4.5f}".format( num_atom, num_mol, coords[i,0], coords[i,1], coords[i,2] )
        num_atom = num_atom + 1
        if ((i+1)%6 == 0):
            num_mol = num_mol + 1

    # start defining how atoms are connected
    print "\nBonds\n"

    # define C-C and C=C bonds here
    for i in range(n_molecules):
        for j in range(n_bonds):
            first = i * n_chain + j + 1
            print "{}	1	{}	{}".format( i * n_bonds + j + 1, first, first + 1 )
  
    # define Cd-CH2 bonds here
    for i in range(n_molecules):
        carbon = i * n_chain + 1
        cadmium = n_molecules * n_chain + i * 2 + 1
        print "{}	1	{}	{}".format( n_molecules * n_bonds + i + 1, carbon, cadmium )

    num_bond = n_molecules * (n_bonds + 1) + 1

    # define hex bonds here
    for i in range(n_hexmolecules):
        for j in range(n_hexbonds):
            first = n_molecules * n_chain + n_molecules * 16 + i * n_hex + j + 1
            print "{}	1	{}	{}".format( num_bond, first, first + 1 )
            num_bond = num_bond + 1

    nat = n_molecules * n_chain + n_molecules * 16 + n_hexatoms + 1
    # define bottom ligand bonds here
    for i in range(n_bottommolecules):
        for j in range(n_bonds):
            print "{}	1	{}	{}".format( num_bond, nat, nat + 1 )
            num_bond = num_bond + 1
            nat = nat + 1
            if (j == n_bonds - 1):
                nat = nat + 1

    # define bottom hexane bonds here
    for i in range(n_hexmolecules):
        for j in range(n_hexbonds):
            print "{}	1	{}	{}".format( num_bond, nat, nat + 1 )
            num_bond = num_bond + 1
            nat = nat + 1
            if (j == n_hexbonds - 1):
                nat = nat + 1


    print "\nAngles\n"

    num_angle = 1
    # oleic acid angles here
    for i in range(n_molecules):
        for j in range(n_angles):
            first = i * n_chain + j + 1
            print "{}	1	{}	{}	{}".format( num_angle, first, first + 1, first + 2 )
            num_angle = num_angle + 1

    # hex angles defined here
    for i in range(n_hexmolecules):
        for j in range(n_hexangles):
            first = n_molecules * n_chain + n_molecules * 16 + i * n_hex + j + 1
            print "{}	1	{}	{}	{}".format( num_angle, first, first + 1, first - 3)
            num_angle = num_angle + 1

    # bottom angles here
    nat = n_molecules * n_chain + n_molecules * 16 + n_hexatoms + 1
    for i in range(n_bottommolecules):
        for j in range(n_angles):
            print "{}	1	{}	{}	{}".format( num_angle, nat, nat + 1, nat + 2 )
            nat = nat + 1
            num_angle = num_angle + 1
            if (j == n_angles - 1):
                nat = nat + 2
    
    # bottom hexane angles here
    for i in range(n_hexmolecules):
        for j in range(n_hexangles):
            print "{}	1	{}	{}	{}".format( num_angle, nat, nat + 1, nat + 2 )
            nat = nat + 1
            num_angle = num_angle + 1
            if (j == n_hexangles - 1):
                nat = nat + 2

    print "\nDihedrals\n"

    num_dihedral = 1
    # oleic acid dihedrals here
    for i in range(n_molecules):
        for j in range(n_dihedrals):
            first = i * n_chain + j + 1
            print "{}	1	{}	{}	{}	{}".format( num_dihedral, first, first + 1, first + 2, first + 3 )
            num_dihedral = num_dihedral + 1

    # hex dihedrals here
    for i in range(n_hexmolecules):
        for j in range(n_hexdihedrals):
            first = n_molecules * n_chain + n_molecules * 16 + i * n_hex + j + 1
            print "{}	1	{}	{}	{}	{}".format( num_dihedral, first, first + 1, first + 2, first + 3 )
            num_dihedral = num_dihedral + 1

    # bottom dihedrals here
    nat = n_molecules * n_chain + n_molecules * 16 + n_hexatoms + 1
    for i in range(n_bottommolecules):
        for j in range(n_dihedrals):
            print "{}	1	{}	{}	{}	{}".format( num_dihedral, nat, nat + 1, nat + 2, nat + 3 )
            num_dihedral = num_dihedral + 1
            nat = nat + 1
            if (j == n_dihedrals - 1):
                nat = nat + 3

    # bottom hexane dihedrals here
    for i in range(n_hexmolecules):
        for j in range(n_hexdihedrals):
            print "{}	1	{}	{}	{}	{}".format( num_dihedral, nat, nat + 1, nat + 2, nat + 3 )
            num_dihedral = num_dihedral + 1
            nat = nat + 1
            if (j == n_hexdihedrals - 1):
                nat = nat + 3


