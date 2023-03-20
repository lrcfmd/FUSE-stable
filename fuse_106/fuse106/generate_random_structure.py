import os
import sys
from fuse106.create_random_string import *
from fuse106.error_check_structure import *
from fuse106.assemble_structure import *
from fuse106.create_random_instructions import *
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
	while complete == False:
			#try:
			#print("\nstart")
			#generate initial string component and start of instructions
			try:
				string,instructions=create_random_string(cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,monoclinic_solutions,atoms_per_fu,fu,vac_ratio=vac_ratio,max_fus=max_fus,system_type=system_type,composition=composition,ap=ap)
			except:
				continue
			#check to see if it has the correct number of atoms
			#print("string: \n",string)
			
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

