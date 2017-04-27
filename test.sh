#!/bin/bash --login

cat << LAMMPSIN > lmp_test.in
# Oleic acid simulation

# initialization
units			real
atom_style		molecular
boundary		p p p
bond_style		harmonic
angle_style		harmonic
dihedral_style		opls
improper_style		harmonic

variable		TargetT equal 370

read_data		all.lmp

group			CH3 type 1
group			CH2single type 2
group			CH2double type 3
group			double type 4
group			ring type 5 6 7
group			noring/solv type 1 2 3 4 10 11
group			solvent type 10 11
group			metal type 8 9
group			nometal type 1 2 3 4 5 6 7 10 11

compute			MyTemp nometal temp

# define interactions
pair_style		lj/cut 14.0
pair_modify		shift yes mix arithmetic

pair_coeff		1 1 0.195 3.75			# for CH3 atoms
pair_coeff      	1 8 0.143  3.54
pair_coeff      	1 9 0.143  3.54

pair_coeff		2 2 0.0914 3.95			# for CH2 atoms with single bonds
pair_coeff      	2 8 0.111  3.54    
pair_coeff      	2 9 0.111  3.54    

pair_coeff		3 3 0.0914 3.95			# for CH2 atoms with adjacent double bond atoms
pair_coeff      	3 8 0.143  3.54
pair_coeff      	3 9 0.143  3.54

pair_coeff		4 4 0.0934 3.73			# for CH=CH double bonds
pair_coeff      	4 8 0.143  3.54
pair_coeff      	4 9 0.143  3.54

pair_coeff		5 5 0.042 3.88			# for C atom at ring intro
pair_coeff      	5 8 0.143  3.54
pair_coeff      	5 9 0.143  3.54

pair_coeff		6 6 0.06 3.7			# for C atoms in ring with no hydrogens
pair_coeff      	6 8 0.143  3.54
pair_coeff      	6 9 0.143  3.54

pair_coeff		7 7 0.101 3.695			# for C-H atoms in ring
pair_coeff      	7 8 0.143  3.54
pair_coeff      	7 9 0.143  3.54

pair_coeff      	9 9 0.0    0.0			# Cd
pair_coeff      	9 10 0.143 3.54
pair_coeff      	9 11 0.111 3.54

pair_coeff      	8 8 0   0.0			# S
pair_coeff      	8 10 0.143 3.54
pair_coeff      	8 11 0.111 3.54

pair_coeff      	10 10 0.195 3.75		# hexane
pair_coeff      	11 11 0.0914 3.95		# hexane

bond_coeff		1 95.9 1.54
bond_coeff		2 95.9 1.34
bond_coeff		3 2000.0 1.40			# ring atoms don't move

angle_coeff		1 62.1 114
angle_coeff		2 70.0 119.7
angle_coeff		3 2000.0 120.0			# ring angles don't change

dihedral_coeff		1 1.4114 -0.2711 3.1458 0
dihedral_coeff		2 0.3432 -0.4363 -1.1217 2.7364
dihedral_coeff		3 0 49.3 0.0 0

improper_coeff		1 2000.0 0.0

# special_bonds  		lj 0.0 0.0 0.0

# simulation parameters
variable		xlo equal 30
variable		xhi equal 75
variable		ylo equal 45
variable		yhi equal 70
variable		zlo equal 25
variable		zhi equal 75
region			BulkSolv block \${xlo} \${xhi} \${ylo} \${yhi} \${zlo} \${zhi} side in units box
variable		SolvDens equal count(solvent,BulkSolv)/((v_xhi-v_xlo)*(v_yhi-v_ylo)*(v_zhi-v_zlo))
variable		TargetDens equal 0.0246

thermo_style		custom step temp c_MyTemp etotal ke pe press v_SolvDens
thermo			10

timestep		1
fix			1 nometal nvt temp \${TargetT} \${TargetT} 50
fix_modify		1 temp MyTemp
dump			traj all xyz 500 dump.test.xyz

velocity		all create \${TargetT} 87287 loop geom
# minimize		1.0e-7 1.0e-6 1000 1000
run			5000
write_restart		restart.test
undump			traj
unfix			1
LAMMPSIN

mpirun -n 8 lmp_mpi -in lmp_test.in
