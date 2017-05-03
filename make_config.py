import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-nx", type=int, default=20, help="number of ligands in x-direction")
parser.add_argument("-nz", type=int, default=12, help="number of ligands in z-direction")
parser.add_argument("-fmt", type=str, choices=["xyz", "lmp"], default='xyz', help="output format")
parser.add_argument("-nolig", action='store_true', help="to make configuration with just dye on the surface")
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
n_dye = 24
n_dyebonds = 27
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
    if args.nolig:
        print "{}".format( n_molecules * 16 + 2 * n_hexatoms + n_dye )
    else:
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

    if not (args.nolig):
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

    if not (args.nolig):
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
    f = open("surface-CdS.txt","w")

    # beginning of LAMMPS input file
    print "LAMMPS coords for {} oleic acid-like molecules\n".format( n_molecules )
    if args.nolig:
        print "{}	atoms".format( n_molecules * 16 + 2 * n_hexatoms + n_dye )
        print "{}	bonds".format( 2 * n_hexmolecules * n_hexbonds + n_dyebonds )
        print "{}	angles".format( 2 * n_hexmolecules * n_hexangles + n_dyeangles )
        print "{}	dihedrals".format( 2 * n_hexmolecules * n_hexdihedrals + n_dyedihedrals )
        print "{} impropers".format( n_dyeimpropers )

    else:
        print "{}	atoms".format( (n_molecules - 5) * n_chain + n_molecules * 16 + 2 * n_hexatoms + n_bottom + n_dye )
        print "{}	bonds".format( (n_molecules - 5) * n_bonds + (n_molecules - 5) + n_molecules / 2 + 2 * n_hexmolecules * n_hexbonds + n_bottommolecules * n_bonds + n_dyebonds )
        print "{}	angles".format( (n_molecules - 5) * n_angles + 2 * n_hexmolecules * n_hexangles + n_bottommolecules * n_angles + n_dyeangles )
        print "{}	dihedrals".format( (n_molecules - 5) * n_dihedrals + 2 * n_hexmolecules * n_hexdihedrals + n_bottommolecules * n_dihedrals + n_dyedihedrals )
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
    print "3 14.027"
    print "4 13.019"
    print "5 12.011"
    print "6 12.011"
    print "7 13.019"
    print "8 32.065"
    print "9 112.411"
    print "10 15.034"
    print "11 14.026\n"
    print "Atoms\n"

    nmol = 1
    natom = 1
    # begin defining actual coordinates

    # dye on surface
    for i in range(n_dye):
        if (i == 6):
            print "{}	{}	3	{:4.5f}	{:4.5f}	{:4.5f}".format( natom, nmol, centred_dye[i,0], centred_dye[i,1], centred_dye[1,2] )
        elif (i == 7): 
            print "{}	{}	4	{:4.5f}	{:4.5f}	{:4.5f}".format( natom, nmol, centred_dye[i,0], centred_dye[i,1], centred_dye[1,2] )
        elif (i == 8):
            print "{}	{}	5	{:4.5f}	{:4.5f}	{:4.5f}".format( natom, nmol, centred_dye[i,0], centred_dye[i,1], centred_dye[1,2] )
        elif (i == 12 or i == 13 or i ==9 or i == 10):
            print "{}	{}	6	{:4.5f}	{:4.5f}	{:4.5f}".format( natom, nmol, centred_dye[i,0], centred_dye[i,1], centred_dye[1,2] )
        elif (i > 8 and i < 22):
            print "{}	{}	7	{:4.5f}	{:4.5f}	{:4.5f}".format( natom, nmol, centred_dye[i,0], centred_dye[i,1], centred_dye[1,2] )
        else:
            print "{}	{}	2	{:4.5f}	{:4.5f}	{:4.5f}".format( natom, nmol, centred_dye[i,0], centred_dye[i,1], centred_dye[1,2] )
        natom = natom + 1

    if not (args.nolig):
        for i in range(nx):
            for j in range(nz):
                coords = np.copy(centred_data)
                # shift molecule
                coords[:, 0] = coords[:, 0] + 4.12 * i
                coords[:, 2] = coords[:, 2] + 6.75 * j

                if (i >= 10 and i <= 14 and j == 6):
                    continue
            
                for k in range(n_chain):
                    if (k == n_chain-1):
                        print "{}	{}	1	{:4.5f}	{:4.5f}	{:4.5f}".format( natom, nmol+1, coords[k,0], coords[k,1], coords[k,2]  )
                    else:
                        print "{}	{}	2	{:4.5f}	{:4.5f}	{:4.5f}".format( natom, nmol+1, coords[k,0], coords[k,1], coords[k,2]  )
                    natom = natom + 1
                nmol = nmol + 1

    # Cd-S surface coordinates
    for i in range(nx):
        for j in range(nz):
            coords = centred_data[0, :] - np.array([0.0, 1.56, 0.0])
            coords[0] = coords[0] + 4.12 * i
            coords[2] = coords[2] + 6.75 * j
            print "{}	{}	9	{:4.5f}	{:4.5f}	{:4.5f}".format( natom, nmol+1, coords[0], coords[1], coords[2] )
            print "{}	{}	8	{:4.5f}	{:4.5f}	{:4.5f}".format( natom + 1, nmol+1, coords[0], coords[1], coords[2] + 2.5299 )
            if (i >= 10 and i <= 14 and j == 6):
                f.write("Cd {} {} {}\n".format(natom, i, j))
            natom = natom + 2

    for i in range(nx):
        for j in range(nz):
            coords = centred_data[0, :] - np.array([0.0, 1.56, 0.0])
            coords[0] = coords[0] + 4.12 * i
            coords[2] = coords[2] + 6.75 * j
            print "{}	{}	9	{:4.5f}	{:4.5f}	{:4.5f}".format( natom,  nmol+1, coords[0] - 2.06073, coords[1] - 1.18977, coords[2] + 3.375 )
            print "{}	{}	8	{:4.5f}	{:4.5f}	{:4.5f}".format( natom + 1, nmol+1, coords[0] - 2.06073, coords[1] - 1.18977, coords[2] - 3.375 + 2.5299 )
            print "{}	{}	9	{:4.5f}	{:4.5f}	{:4.5f}".format( natom + 2, nmol+1, coords[0] - 2.06073, coords[1] - 3.5693, coords[2] )
            print "{}	{}	8	{:4.5f}	{:4.5f}	{:4.5f}".format( natom + 3, nmol+1, coords[0] - 2.06073, coords[1] - 3.5693, coords[2] + 2.5299 )
            print "{}	{}	9	{:4.5f}	{:4.5f}	{:4.5f}".format( natom + 4, nmol+1, coords[0], coords[1] - 3.5693 - 1.18977, coords[2] + 3.375 )
            print "{}	{}	8	{:4.5f}	{:4.5f}	{:4.5f}".format( natom + 5, nmol+1, coords[0], coords[1] - 3.5693 - 1.18977, coords[2] - 3.375 + 2.5299 )

            print "{}	{}	9	{:4.5f}	{:4.5f}	{:4.5f}".format( natom + 6, nmol+1, coords[0], coords[1] - 7.13848, coords[2] )
            print "{}	{}	8	{:4.5f}	{:4.5f}	{:4.5f}".format( natom + 7, nmol+1, coords[0], coords[1] - 7.13848, coords[2] + 2.5299 )
            print "{}	{}	9	{:4.5f}	{:4.5f}	{:4.5f}".format( natom + 8, nmol+1, coords[0] - 2.06073, coords[1] - 7.13848 - 1.18977, coords[2] + 3.375 )
            print "{}	{}	8	{:4.5f}	{:4.5f}	{:4.5f}".format( natom + 9, nmol+1, coords[0] - 2.06073, coords[1] - 7.13848 - 1.18977, coords[2] - 3.375 + 2.5299 )
            print "{}	{}	9	{:4.5f}	{:4.5f}	{:4.5f}".format( natom + 10, nmol+1, coords[0] - 2.06073, coords[1] - 7.13848 - 3.5693, coords[2] )
            print "{}	{}	8	{:4.5f}	{:4.5f}	{:4.5f}".format( natom + 11, nmol+1, coords[0] - 2.06073, coords[1] - 7.13848 - 3.5693, coords[2] + 2.5299 )
            print "{}	{}	9	{:4.5f}	{:4.5f}	{:4.5f}".format( natom + 12, nmol+1, coords[0], coords[1] - 7.13848 - 3.5693 - 1.18977, coords[2] + 3.375 )
            print "{}	{}	8	{:4.5f}	{:4.5f}	{:4.5f}".format( natom + 13, nmol+1, coords[0], coords[1] - 7.13848 - 3.5693 - 1.18977, coords[2] - 3.375 + 2.5299 )
            natom = natom + 14

    # hex molecules here

    nmol = nmol + 2
    coords = centred_hex + np.array([0.0, 0.0, 0.0]) 
    for i in range(n_hexatoms):
        if ((i%6 == 0) or ((i+1)%6 == 0)):
            print "{}	{}	10	{:4.5f}	{:4.5f}	{:4.5f}".format( natom, nmol, coords[i,0], coords[i,1], coords[i,2] )
        else:
            print "{}	{}	11	{:4.5f}	{:4.5f}	{:4.5f}".format( natom, nmol, coords[i,0], coords[i,1], coords[i,2] )

        natom = natom + 1
        # increase molecule number every 6 atoms
        if ((i+1)%6 == 0):
            nmol = nmol + 1

    if not (args.nolig):
        # bottom molecules here
        coords = centred_bottom + np.array([0.0, 0.0, 0.0])
        for i in range(n_bottom):
            if ((i % (n_chain-1)) == 0 and i != 0):
                print "{}	{}	1       {:4.5f} {:4.5f} {:4.5f}".format( natom, nmol, coords[i,0], coords[i,1] + -15.457550, coords[i,2] )
            else:
                print "{}	{}	2       {:4.5f} {:4.5f} {:4.5f}".format( natom, nmol, coords[i,0], coords[i,1] + -15.457550, coords[i,2] )
            natom = natom + 1
            if (i%n_chain == 0):
                nmol = nmol + 1

    # bottom hex molecules here
    coords = centred_hex + np.array([0.0, -55.5 - 80, 0.0])
    for i in range(n_hexatoms):    
        if ((i%6 == 0) or ((i+1)%6 == 0)):
            print "{}	{}	10	{:4.5f}	{:4.5f}	{:4.5f}".format( natom, nmol, coords[i,0], coords[i,1], coords[i,2] )
        else:
            print "{}	{}	11	{:4.5f}	{:4.5f}	{:4.5f}".format( natom, nmol, coords[i,0], coords[i,1], coords[i,2] )
        natom = natom + 1
        if ((i+1)%6 == 0):
            nmol = nmol + 1

    f.close()
    # start defining how atoms are connected
    print "\nBonds\n"

    # define dye bonds here
    print "1       1       1       2"
    print "2       1       2       3"
    print "3       1       3       4"
    print "4       1       4       5"
    print "5       1       5       6"
    print "6       1       6       7"
    print "7       1       7       8"
    print "8       2       8       9"
    print "9       3       9       10"
    print "10      3       10      11" 
    print "11      3       11      12" 
    print "12      3       12      13" 
    print "13      3       13      14" 
    print "14      3       14      9" 
    print "15      3       11      16"  
    print "16      3       16      17" 
    print "17      3       17      18" 
    print "18      3       18      19" 
    print "19      3       19      10" 
    print "20      3       13      15" 
    print "21      3       15      20" 
    print "22      3       20      21"
    print "23      3       21      22" 
    print "24      3       22      14" 
    print "25      1       1       23"
    print "26      1       23      24" 

    nbond = 28

    if not (args.nolig):
        print "27      1       24      4507"
        # define ligand bonds here
        for i in range(n_molecules - 5):
            for j in range(n_bonds):
                first = i * n_chain + j + 1 + n_dye
                print "{}	1	{}	{}".format( nbond, first, first + 1 )
                nbond = nbond + 1
      
        # define top Cd-CH2 bonds here
        cadmium = (n_molecules - 5) * n_chain + 1 + n_dye
        for i in range(n_molecules - 5):
            carbon = i * n_chain + 1 + n_dye
            # check if ligand was removed
            if (cadmium == 4507 or cadmium == 4531 or cadmium == 4555 or cadmium == 4579 or cadmium == 4603):
                cadmium = cadmium + 2
            print "{}	1	{}	{}".format( nbond, carbon, cadmium )
            cadmium = cadmium + 2
            nbond = nbond + 1
    
        nat = (n_molecules - 5) * n_chain + n_dye + n_molecules * 16 + 1
    else:
        print "27      1       24      277"
        nat = n_molecules * 16 + 1

    # define hex bonds here
    for i in range(n_hexmolecules):
        for j in range(n_hexbonds):
            print "{}	1	{}	{}".format( nbond, nat, nat + 1 )
            nbond = nbond + 1
            nat = nat + 1
            if (j == n_hexbonds - 1):
                nat = nat + 1

    if not (args.nolig):
        # define bottom Cd-CH2 bonds here
        cadmium = cadmium + 12
        lig_in_row = 0
        for i in range(n_bottommolecules):
            carbon = nat + n_chain * i
            print "{}	1	{}	{}".format( nbond, carbon, cadmium )
            cadmium = cadmium + 14
            lig_in_row = lig_in_row + 1
            nbond = nbond + 1
            # skip full row of cadmium atoms w/o ligands
            if (lig_in_row == 12):
                cadmium = cadmium + 12 * 14
                lig_in_row = 0

       # define bottom ligand bonds here
        for i in range(n_bottommolecules):
            for j in range(n_bonds):
                print "{}	1	{}	{}".format( nbond, nat, nat + 1 )
                nbond = nbond + 1
                nat = nat + 1
                if (j == n_bonds - 1):
                    nat = nat + 1

    # define bottom hexane bonds here
    for i in range(n_hexmolecules):
        for j in range(n_hexbonds):
            print "{}	1	{}	{}".format( nbond, nat, nat + 1 )
            nbond = nbond + 1
            nat = nat + 1
            if (j == n_hexbonds - 1):
                nat = nat + 1

    print "\nAngles\n"

    # print dye bonds here
    print "1       1       1       2       3"
    print "2       1       2       3       4"
    print "3       1       3       4       5"
    print "4       1       4       5       6"
    print "5       1       6       7       8"
    print "6       2       7       8       9"
    print "7       3       8       9       10"
    print "8       3       8       9       14"
    print "9       3       9       10      11" 
    print "10      3       10      11      12"
    print "11      3       11      12      13"
    print "12      3       12      13      14" 
    print "13      3       13      14      9"
    print "14      3       10      11      16"
    print "15      3       11      16      17" 
    print "16      3       16      17      18"
    print "17      3       17      18      19" 
    print "18      3       18      19      10"  
    print "19      3       14      13      15" 
    print "20      3       13      15      20" 
    print "21      3       15      20      21" 
    print "22      3       20      21      22"
    print "23      3       21      22      14" 
    print "24      3       22      14      9" 
    print "25      3       14      9       10"
    print "26      3       9       10      19" 
    print "27      3       15      13      12"
    print "28      3       12      11      16"
    print "29      1       24      23      1"
    print "30      1       23      1       2"

    num_angle = 31

    if not (args.nolig):
        # ligand angles here
        for i in range(n_molecules - 5):
            for j in range(n_angles):
                first = i * n_chain + j + 1 + n_dye
                print "{}	1	{}	{}	{}".format( num_angle, first, first + 1, first + 2 )
                num_angle = num_angle + 1
    
        nat = (n_molecules - 5) * n_chain + n_dye + n_molecules * 16 + 1
    else:
        nat = n_molecules * 16 + 1

    # hex angles defined here
    for i in range(n_hexmolecules):
        for j in range(n_hexangles):
            print "{}	1	{}	{}	{}".format( num_angle, nat, nat + 1, nat + 2)
            num_angle = num_angle + 1
            nat = nat + 1
            if (j == n_angles - 1):
                nat = nat + 2

    if not (args.nolig):
        # bottom angles here
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

    # dye dihedrals here
    print "1       1       1       2       3       4"
    print "2       1       2       3       4       5"
    print "3       1       3       4       5       6"
    print "4       1       4       5       6       7"
    print "5       1       5       6       7       8"
    print "6       2       6       7       8       9"
    print "7       1       24      23      1       2"
    print "8       1       23      1       2       3"
    print "9       3       7       8       9       10"

    num_dihedral = 10

    if not (args.nolig):
        # ligand dihedrals here
        for i in range(n_molecules - 5):
            for j in range(n_dihedrals):
                first = i * n_chain + j + 1 + n_dye
                print "{}	1	{}	{}	{}	{}".format( num_dihedral, first, first + 1, first + 2, first + 3 )
                num_dihedral = num_dihedral + 1
    
        nat = (n_molecules - 5) * n_chain + n_dye + n_molecules * 16 + 1
    else:
        nat = n_molecules * 16 + 1

    # hex dihedrals here
    for i in range(n_hexmolecules):
        for j in range(n_hexdihedrals):
            print "{}	1	{}	{}	{}	{}".format( num_dihedral, nat, nat + 1, nat + 2, nat + 3 )
            num_dihedral = num_dihedral + 1
            nat = nat + 1
            if (j == n_hexdihedrals - 1):
                nat = nat + 3

    if not (args.nolig):
        # bottom dihedrals here
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

    print "\nImpropers\n"

    # dye impropers here
    print "1       1       9       8       10      14"
    print "2       1       14      9       22      13"
    print "3       1       13      14      15      12"
    print "4       1       10      9       11      19"
    print "5       1       11      12      10      16\n"

