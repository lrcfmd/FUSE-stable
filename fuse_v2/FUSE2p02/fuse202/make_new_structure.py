#fuse imports
from fuse202.possible_solutions import *
from fuse202.generate_random_structure import *
from fuse202.extract_modules import *
import fuse202.bond_table

#other imports
from random import choice
import math
import numpy
from ase.io import *
from ase.visualize import *
import pandas

def get_new_structure(composition='',max_atoms='',imax_atoms='',restart='',
	max_ax='',density_cutoff ='',check_bonds='',btol='',check_distances='',
	system_type='', dist_cutoff = '',vac_ratio = '',atoms_per_fu='',imax_fus='',
	max_fus='',cubic_solutions='',tetragonal_solutions='',hexagonal_solutions='',
	orthorhombic_solutions='',monoclinic_solutions='',bondtable='',
	ideal_density='',fu='',ap='',use_spglib=''
):
	gen = False
	while gen == False:
		##############################################################################
		# now call the generate random structure function
		#for each iteration of generating n structures:
		target_fu=choice(list(range(1,imax_fus+1)))
		target_atoms=target_fu*atoms_per_fu
		#print("hello")
		#print(target_fu)
		#print(target_atoms)
		atoms,string,instructions,accept=generate_random_structure(target_atoms,cubic_solutions,tetragonal_solutions, hexagonal_solutions, orthorhombic_solutions, monoclinic_solutions, atoms_per_fu, ideal_density, density_cutoff, check_bonds, btol, system_type, fu, composition, bondtable, ap, check_distances, dist_cutoff, vac_ratio=vac_ratio, max_fus=imax_fus, target_number_atoms=target_atoms)
		
		##############################################################################
		# now call write the atoms object to a temporary file & then call the extract modules function
		if use_spglib == True:
			import spglib
			from ase import Atoms
			try:
				lattice,positions,numbers=spglib.standardize_cell(atoms,symprec=1.e-5)
				temp2=Atoms(numbers=numbers,pbc=True)
				temp2.cell=lattice
				temp2.set_scaled_positions(positions)
				temp2.numbers=numbers
				atoms=temp2.copy()
				
				#print(atoms.cell.cellpar())
				
			except:
				pass
   	
		if not atoms == None:
			write("temp.cif",atoms)
			input_files=["temp.cif"]
			structure=extract_module(input_files,bondtable)
			gen = True

	#print(structure['sub module cell'])
	#print(structure['shape in submods'])
	
	return structure