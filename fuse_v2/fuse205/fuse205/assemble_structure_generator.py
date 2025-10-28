from ase import *
from ase.visualize import *
from random import shuffle
import sys
import warnings
import os
# function to create atoms object from a set of assembly instructions and a string
# of atomic numbers
def assemble_structure(string,instructions):
	temp_string=string.copy()
	### reading assembly instructions ###########################################
	ap=instructions[0] # lattce parameter
	lattice=instructions[1] # lattice type: 
	full=instructions[2:5] # unit cell dimensions in units of sub-modules
	angles=instructions[5:8] # unit cell angles
	translations=instructions[8:] #translation instructions
	# 0 = cubic
	# 1 = tetragonal
	# 2 = hexagonal
	# 3 = orthorhombic
	# 4 = monoclinic
	# 5 = triclinic
	#if os.path.isfile("lattice.txt"):
	#	f=open("lattice.txt",'a')
	#	f.write(str(lattice)+str("\n"))
	#	f.close()
	#else:
	#	f=open("lattice.txt",'w')
	#	f.write(str(lattice)+str("\n"))
	#	f.close()
	#############################################################################
	
	### splice the input string into sub-modules	################################
	nsub=int(len(temp_string)/4)
	sub_mods=[]
	for i in range(nsub):
		temp=temp_string[-4:]
		sub_mods.append(temp)
		del temp_string[-4:]
	#############################################################################
	
	### make the submodules into atoms objects ##################################
	co_ords=[[0,0,0.5],[0.5,0,0.5],[0.5,0.5,0.5],[0,0.5,0.5]]
	sub_mods_atoms=[]
	for i in range(nsub):
		atoms=None
		for j in range(4):
			if sub_mods[i][j] !=120:
				tmp_atoms=Atoms(numbers=[sub_mods[i][j]],pbc=[1,1,1],cell=[ap,ap,ap*0.5,angles[0],angles[1],angles[2]],scaled_positions=[co_ords[j]])
				if atoms != None:
					atoms = atoms + tmp_atoms
				if atoms == None:
					atoms=tmp_atoms
					
		if atoms == None:
			atoms=Atoms(pbc=[1,1,1],cell=[ap,ap,0.5*ap,angles[0],angles[1],angles[2]])
		sub_mods_atoms.append(atoms)
		
	
	#############################################################################
	
	### assemble the sub-modules into modules ###################################
	mod_size=[full[0],full[1]]
	mods=[]
	nsub_mods=full[0]*full[1]
	for x in range(full[2]):
		n=0
		temp=sub_mods_atoms[-nsub_mods:]
		for i in range(mod_size[0]):
			for j in range(mod_size[1]):
				if i + j == 0:
					atoms=temp[0]
					n+=1
				else:
					temp_atoms=temp[n]
					n+=1
					temp_atoms.translate(atoms.cell[0]*i)
					temp_atoms.translate(atoms.cell[1]*j)
					atoms+=temp_atoms
				
		atoms.cell[0]=atoms.cell[0]*full[0]
		atoms.cell[1]=atoms.cell[1]*full[1]
		del sub_mods_atoms[-nsub_mods:]
		mods.append(atoms)
	
	#############################################################################

	# define which translations are in use ######################################
	if angles[-1] == 90:
		ptrans=[[0.,0.,0.],[0.5,0.5,0.]]
	if angles[-1] == 120:
		ptrans=[[0,0,0],[1/3.,2/3.,0],[2/3.,1/3.,0]]		
	### then need to stich the modules together into a full structure ###########
	for i in range(full[2]):
		if i == 0:
			atoms = mods[i]
			atoms.translate(([ptrans[translations[i]][0]*ap,ptrans[translations[i]][1]*ap,ptrans[translations[i]][2]*(ap/2)]))
		else:
			temp_atoms=mods[i]
			temp_atoms.translate(atoms.cell[2])
			temp_atoms.translate(([ptrans[translations[i]][0]*ap,ptrans[translations[i]][1]*ap,ptrans[translations[i]][2]*(ap/2)]))
			atoms+=temp_atoms
			atoms.cell[2]=atoms.cell[2]+temp_atoms.cell[2]
	atoms.set_pbc=([True,True,True])
	
	#############################################################################
	
	# now need to set the unit cell angles
	#atoms.cell=[atoms.get_cell_lengths_and_angles()[0],atoms.get_cell_lengths_and_angles()[1],atoms.get_cell_lengths_and_angles()[2],angles[0],angles[1],angles[2]]
	#write("atoms.cif",atoms)
	return atoms, instructions
	
