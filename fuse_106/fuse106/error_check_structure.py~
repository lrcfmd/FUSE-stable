from fuse106.get_distances import *
import sys
from fuse106.test_bonds import test_bonds
import numpy
from ase.visualize import *

def error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,
	system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,check_anions=False,
	target_number_atoms=''):
	accept=1 # accept the structure unless told otherwise
	###
	#write("test.cif",atoms)
	###
	### check to see if atoms object has physical unit cell ########
	cell=atoms.get_cell_lengths_and_angles()
	if numpy.isnan(cell).any():
		accept=0
	#############################################################################

	#############################################################################
	### check to see if we are within the limit on the number of atoms 
	n_at=len(atoms)
	if n_at > target_number_atoms:
		accept = 0
		return accept
	
	### now need to explicitly check the composition agaisnt the input composition
	#print(composition)
	#sys.exit()

	### check unit cell density #################################################
	mass=atoms.get_masses().sum()
	volume=atoms.get_volume()
	density=mass/volume
	#ideal_density = float(ideal_density)
	#density_cutoff = float(density_cutoff)
	#print(density_cutoff)
	cut=float(ideal_density) * float(density_cutoff)

	if density < cut:
		accept=0
		return accept
	if density > ideal_density:
		accept=0		
		return accept
	
	
	#############################################################################
	
	### check for the shortest distance in the structure, if any less than cutoff
	# distance, reject structure 
	if check_distances==True:
		temp_atoms=atoms.repeat([2,2,2])
		#distances=min(get_distances(new_atoms=atoms))
		#print (distances)
		temp1=temp_atoms.get_all_distances()
		temp2=[]
		for i in range(len(temp1)):
			for j in range(len(temp1[i])):
				if temp1[i][j] != 0:
					temp2.append(temp1[i][j])
		distances=min(temp2)
		#print(distances,end=' ')
		if distances <= dist_cutoff:
			accept=0
			return accept
			
	#############################################################################
	
	### if in use, examine the bonding in the proposed structure ################
	# If neutral input, will have to test all co-ordinations, if ionic, just use 
	# the cation species as the centres
	if check_bonds == True:
		if system_type=='neutral':
			cations=[]
			anions=[]
			charges=[]
			for i in range(len(fu)):
				temp=Atoms(numbers=[fu[i]]).get_chemical_symbols()[0]
				temp2=composition[temp]
				if not temp in anions:
					anions.append(temp)
				if not temp in cations:
					cations.append(temp)
			
			errors=test_bonds(atoms=atoms,cations=cations,anions=anions,charges=charges,ap=ap,lib=bondtable,system_type=system_type)
			if errors > btol:
				accept=0
				return accept
				
		if system_type=='ionic':
			cations=[]
			anions=[]
			charges=[]
			an_charges=[]
			for i in range(len(fu)):
				temp=Atoms(numbers=[fu[i]]).get_chemical_symbols()[0]
				temp2=composition[temp]
				if temp2[1] < 0:
					if not temp in anions:
						anions.append(temp)
				if temp2[1] > 0:
					if not temp in cations:
						cations.append(temp)
			
			for i in range(len(cations)):
				charges.append(composition[cations[i]][1])
			
			for i in range(len(anions)):
				an_charges.append(composition[anions[i]][1])
			
			errors=test_bonds(atoms=atoms,cations=cations,anions=anions,charges=charges,ap=ap,lib=bondtable,system_type=system_type)
			#print("errors 1: " + str(errors))
			
			if check_anions==True:
				an_errors=test_bonds(atoms=atoms,cations=anions,anions=cations,charges=an_charges,ap=ap,lib=bondtable,system_type=system_type)
				#print("errors 2: " + str(an_errors))
				total_errors=0
				cats=0
				anis=0
				for i in range(len(cations)):
					cats+=atoms.get_chemical_symbols().count(cations[i])
				for i in range(len(anions)):
					anis+=atoms.get_chemical_symbols().count(anions[i])
				cats=cats*errors
				anis=anis*an_errors
				total_errors=cats+anis
				total_errors=total_errors/len(atoms)
				errors=total_errors
			#sys.exit()	
			if errors > btol:
				accept=0
				return accept
			
	#############################################################################
	
	return accept
