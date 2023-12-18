import os
import sys
from fuse202.create_random_string import *
from fuse202.error_check_structure import *
from fuse202.assemble_structure_generator import *
from fuse202.create_random_instructions import *
from ase.visualize import *

def generate_random_structure(target_atoms,cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,
					monoclinic_solutions,atoms_per_fu,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,
					vac_ratio='',max_fus='',
					target_number_atoms=''):

	accept=0
	n_atoms=0
	complete = False
	max_attempts=1000
	attempts = 0
	
	#print(cubic_solutions)
	#sys.exit()
	
	while complete == False:
			#print(attempts)
			#try:
			#print("\nstart")
			#generate initial string component and start of instructions
			#before we try this, print out a list of the variables that we're passing to the random string
			#print("cubic_solutions: ",cubic_solutions)
			#print("tetragonal solitions: ",tetragonal_solutions)
			#print("hexagonal_solutions: ",hexagonal_solutions)
			#print("orthorhombic_solutions: ",orthorhombic_solutions)
			#print("monoclinic_solutions: ",monoclinic_solutions)
			#print("atoms_per_fu: ",atoms_per_fu)
			#print("fu: ",fu)
			#print("vac_ratio: ",vac_ratio)
			#print("max_fus: ", max_fus)
			#print("system_type: ", system_type)
			#print("composition: ",composition)
			#print("ap: ",ap)
			
			try:
				string,instructions=create_random_string(cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,monoclinic_solutions,atoms_per_fu,fu,vac_ratio=vac_ratio,max_fus=max_fus,system_type=system_type,composition=composition,ap=ap)
			except:
				continue
			#check to see if it has the correct number of atoms
			#print("string: \n",string)
			#print("hello2")
			instructions=create_random_instructions(string,ap,cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,monoclinic_solutions,instructions)
			#print("instructions: \n",instructions)
			
			#print("atoms: \n",atoms)
			
			try:
				
				stripped=list(filter((120).__ne__,string))
				n_atoms=len(stripped)
				
			except:
				n_atoms=0
			
			try:
				atoms,instructions=assemble_structure(string,instructions)
			except:
				continue
			
			if len(atoms) > 0:
				accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=target_atoms)
			else:
				accept=0
			
			#print(accept)
			#accept=1
			#except:
			#	pass
			
			if n_atoms == target_atoms:
				if accept == 1:
					complete=True
			
			attempts +=1
			
			if attempts == max_attempts:
				accept = 0
				atoms=None
				string=None
				instructions=None
				break
	#view(atoms)
	#print(n_atoms)
	#print("\n")
	#print("string:\n",string,"\ninstructions: \n",instructions)
	#sys.exit()
	return atoms,string,instructions,accept

