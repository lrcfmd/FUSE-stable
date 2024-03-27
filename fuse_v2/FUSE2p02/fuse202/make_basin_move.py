from ase.atoms import *
from ase.visualize import *
from ase.io import *
from random import choice
from random import shuffle
from fuse202.extract_modules import *
from fuse202.assemble_structure_2 import *
from fuse202.make_new_structure import *
import os
import sys
import math
import random

def make_basin_move(current_structure,moves,bondtable,grid_spacing,exclusion,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,check_distances,dist_cutoff,
	cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,monoclinic_solutions,max_atoms,max_ax,
	vac_ratio,using_prebuilt,pre_built_structures,max_fus,atoms_per_fu,
	imax_atoms,use_spglib,initial_population):
	#set this to keep attempting moves until we make a valid move
	complete = False
	
	#setup the new structure object to return to the main code:
	#set the empty structure dictionary that we're going to populate
	structure={
	'modules':'the module set derrived for the given structure',
	'sub module cell':'shape in [x,y,z,al,be,ga]',
	'shape in submods':'as for sub module cell, but x,y,z are the number of sub modules',
	'nmods':'integer giving the total number of sub modules in the structure',
	'ap':'the root lattice parameter applied to the structure',
	'atoms':'the fully assembled atoms object, this should get updated to the relaxed structure once it has been optimised',
	'optimised?':False,
	'energy':0.0,
	'converged':False
	}
	
	
	trial=0
	
	#for testing... just have a quick look at the current structure
	#view(current_structure['atoms'])
	#write("start.cif",current_structure['atoms'])
	
	# Choose the move that we're going to attempt. I've put it outside of the 
	# while loop so that we should keep attempting the same move until either it works
	# or we meet some condition which forces the move choice to change.
	
	start_structure=current_structure.copy()
	#print("starting atoms: ",len(start_structure['atoms']))
	moves_to_choose=[]
	
	for i in list(moves.keys()):
		for j in range(moves[i]):
			moves_to_choose.append(i)
			
	move=choice(moves_to_choose)
	
	
	
	#*** remember to remove this
	#move = 10
	#***************************
	
	while complete == False:
		print("move: ",move)
		##########################################################################
		if move == 1: #swap the position of two atoms
			atoms=current_structure['atoms'].copy()
			#check that we have more than one element present:
			syms=[]
			for i in atoms:
				if not i.symbol in syms:
					syms.append(i.symbol)
			
			# if we're dealing with an element, redraw a move and continue
			if len(syms) == 1:
				move = move=choice(moves_to_choose)
				continue
			
			if len(syms) > 1:
			
				# First choose two different atoms in the structure & swap their symbol
				
				a1=[choice(list(range(len(atoms))))]
				a2=[choice(list(range(len(atoms))))]			
				a1.append(atoms[a1[0]].symbol)
				a2.append(atoms[a2[0]].symbol)
				
				while a1[1] == a2[1]:
					a2=[choice(list(range(len(atoms))))]			
					a2.append(atoms[a2[0]].symbol)
				
				# now we've got our two atoms, swap them in the structure
				
				atoms[a1[0]].symbol=a2[1]
				atoms[a2[0]].symbol=a1[1]
				
				# before we pass the structure back, check that we have a set of atoms which is consistent with the input composition
				trial_fus=len(atoms)/sum(list(composition.values()))
				counts=[]
				correct=False
				if float(trial_fus).is_integer():
					#now need to check each species
					symbols=atoms.get_chemical_symbols()
					for x in list(composition.keys()):
						num1=symbols.count(x)
						counts.append(num1/composition[x])
					if all(x == trial_fus for x in counts):
						correct = True
				
				if correct == True:
					comp_atoms=atoms.copy()
					complete=True
		
		##########################################################################

		
		##########################################################################
		
		if move == 2: #locate any vacancies in the structure & move an atom into one of them
			# grab a copy of the current atoms object
			atoms=current_structure['atoms'].copy()
			# step one, try and locate vacancies
			# work out what fraction of the cell the grid spacing is in each direction:
			
			
			cell=current_structure['atoms'].cell.cellpar()
			steps=[grid_spacing/cell[0],grid_spacing/cell[1],grid_spacing/cell[2]]
			# initialize our grid
			grid=[] 
			
			#work out how many poitns we need along each axis
			#print(cell)
			npoints=[math.floor(cell[0]/grid_spacing),math.floor(cell[1]/grid_spacing),math.floor(cell[2]/grid_spacing)]
			#populate the grid
			
			for x in range(npoints[0]):
				for y in range(npoints[1]):
					for z in range(npoints[2]):
						grid.append([(x*steps[0])*cell[0],(y*steps[1])*cell[1],(z*steps[2])*cell[2]])
			
			#as a test, make a dummy cell with the grid in:
			full_grid=Atoms("",cell=current_structure['atoms'].cell,pbc=True)
			
			viable=[]
			
			for i in grid:
				at=Atom("X",position=i)
				full_grid+=at
				
				#now work out if the position is viable
				temp_atoms=current_structure['atoms'].copy()
				temp_atoms+=at
				distances=temp_atoms.get_distances(-1,list(range(len(temp_atoms)-1)),mic=True)
				if min(distances) > exclusion:
					viable.append(i)
						
			#now, if there are any viable points, choose a point and an atom and move it
			
			if len(viable) > 0:
				point=choice(viable)
				at=choice(list(range(len(atoms))))
				atoms[at].position=point
				
				# before we pass the structure back, check that we have a set of atoms which is consistent with the input composition
				trial_fus=len(atoms)/sum(list(composition.values()))
				counts=[]
				correct=False
				if float(trial_fus).is_integer():
					#now need to check each species
					symbols=atoms.get_chemical_symbols()
					for x in list(composition.keys()):
						num1=symbols.count(x)
						counts.append(num1/composition[x])
					if all(x == trial_fus for x in counts):
						correct = True
						
				if correct == True:
					comp_atoms=atoms.copy()
					complete=True
				
			#if we have no viable points, redraw the move type and continue:
			if len(viable) == 0:
				move = move=choice(moves_to_choose)
				continue
			
			#view(temp_atoms)
			
		##########################################################################
		
		##########################################################################
			
		if move == 3: # Swap the positions of more than two atoms (but less than all)
			atoms = current_structure['atoms'].copy()
			
			#check to see that there are enough atoms & species to swap in the structure
			
			syms=[]
			for i in atoms:
				if not i.symbol in syms:
					syms.append(i.symbol)
			
			# if we're dealing with an element, redraw a move and continue
			if len(syms) == 1:
				move = choice(moves_to_choose)
				continue			
			
			# we obviously need more than 3 atoms in the structure to make this move
			if len(atoms) <= 3:
				move = choice(moves_to_choose)
				continue			
			
			# then if we have enough atoms / species, get swapping!
			
			#flag to say we're ready to change the atoms objecct
			ready = False
			while ready == False:
				to_swap=choice(list(range(3,len(atoms))))
				
				pool=list(range(len(atoms)))
				chosen=[[choice(pool)]]
				del pool[chosen[0][0]]
				unique_syms=[]
				#print(chosen[0])
				#print(atoms[chosen[0]])
				chosen[0].append(atoms[chosen[0]].symbols[0])
				for i in range(to_swap-1):
					c=[choice(pool)]
					del pool[pool.index(c[0])]
					c.append(atoms[c[0]].symbol)
					chosen.append(c)
					if c[1] not in unique_syms:
						unique_syms.append(c[1])
										
				if len(unique_syms) > 1:
					ready = True
			
			#build a new list with the symbols shuffled
			shuffled=[]
			indicies=[]
			symbols=[]
			
			for i in chosen:
				indicies.append(i[0])
				symbols.append(i[1])
			
			o_symbols=symbols.copy()
			
			while o_symbols == symbols:
				shuffle(symbols)
			
			for i in range(len(indicies)):
				c=[indicies[i],symbols[i]]
				shuffled.append(c)
			
			# now go through and reasign the chemical symbols in the structure			
			for i in shuffled:
				atoms[i[0]].symbol=i[1]
			# before we pass the structure back, check that we have a set of atoms which is consistent with the input composition
			trial_fus=len(atoms)/sum(list(composition.values()))
			counts=[]
			correct=False
			if float(trial_fus).is_integer():
				#now need to check each species
				symbols=atoms.get_chemical_symbols()
				for x in list(composition.keys()):
					num1=symbols.count(x)
					counts.append(num1/composition[x])
				if all(x == trial_fus for x in counts):
					correct = True
			
			if correct == True:
				comp_atoms=atoms.copy()
				complete=True
		##########################################################################
				
		##########################################################################
			
		if move == 5: # Swap the positions of all the atoms in the cell
			atoms = current_structure['atoms'].copy()
			
			#check to see that there are enough atoms & species to swap in the structure
			
			syms=[]
			for i in atoms:
				if not i.symbol in syms:
					syms.append(i.symbol)
			
			# if we're dealing with an element, redraw a move and continue
			if len(syms) == 1:
				move = choice(moves_to_choose)
				continue			
			
			# we obviously need more than 3 atoms in the structure to make this move, if == 3, need to swap all
			if len(atoms) <= 3:
				move = choice(moves_to_choose)
				continue			
			
			# then if we have enough atoms / species, get swapping!
			
			#flag to say we're ready to change the atoms objecct
			ready = False
			while ready == False:
				to_swap=len(atoms)
				
				pool=list(range(len(atoms)))
				chosen=[[choice(pool)]]
				del pool[chosen[0][0]]
				unique_syms=[]
				#print(chosen[0])
				#print(atoms[chosen[0]])
				chosen[0].append(atoms[chosen[0]].symbols[0])
				for i in range(to_swap-1):
					c=[choice(pool)]
					del pool[pool.index(c[0])]
					c.append(atoms[c[0]].symbol)
					chosen.append(c)
					if c[1] not in unique_syms:
						unique_syms.append(c[1])
										
				if len(unique_syms) > 1:
					ready = True
			
			#build a new list with the symbols shuffled
			shuffled=[]
			indicies=[]
			symbols=[]
			
			for i in chosen:
				indicies.append(i[0])
				symbols.append(i[1])
			
			o_symbols=symbols.copy()
			
			while o_symbols == symbols:
				shuffle(symbols)
			
			for i in range(len(indicies)):
				c=[indicies[i],symbols[i]]
				shuffled.append(c)
			
			# now go through and reasign the chemical symbols in the structure
			#write("before.cif",atoms)
			
			for i in shuffled:
				atoms[i[0]].symbol=i[1]
			
			#write("after.cif",atoms)
			#sys.exit()
			# before we pass the structure back, check that we have a set of atoms which is consistent with the input composition
			trial_fus=len(atoms)/sum(list(composition.values()))
			counts=[]
			correct=False
			if float(trial_fus).is_integer():
				#now need to check each species
				symbols=atoms.get_chemical_symbols()
				for x in list(composition.keys()):
					num1=symbols.count(x)
					counts.append(num1/composition[x])
				if all(x == trial_fus for x in counts):
					correct = True
			
			if correct == True:
				comp_atoms=atoms.copy()
				complete=True
				
		##########################################################################
				
		##########################################################################
		# swapping the position of more than two atoms - including vacancies
		if move == 4:	
			# have started by first copying the code from move = 3: swapping some atoms,
			# then need to locate the vacancies within a structure & then decide how many
			# vacancy sites it wants to use when swapping atoms.
			atoms = current_structure['atoms'].copy()
			
			#check to see that there are enough atoms & species to swap in the structure
			
			syms=[]
			for i in atoms:
				if not i.symbol in syms:
					syms.append(i.symbol)
			
			# if we're dealing with an element, redraw a move and continue
			if len(syms) == 1:
				move = choice(moves_to_choose)
				continue			
			
			# we obviously need more than 3 atoms in the structure to make this move, if there are only three, it would need to be the swap all
			if len(atoms) <= 3:
				move = choice(moves_to_choose)
				continue			
			
			#print(len(atoms))
			# then if we have enough atoms / species, get swapping!
			
			#flag to say we're ready to change the atoms objecct
			ready = False
			while ready == False:
				
				to_swap=choice(list(range(3,len(atoms))))
				
				pool=list(range(len(atoms)))
				chosen=[[choice(pool)]]
				del pool[chosen[0][0]]
				unique_syms=[]
				#print(chosen[0])
				#print(atoms[chosen[0]])
				chosen[0].append(atoms[chosen[0]].symbols[0])
				for i in range(to_swap-1):
					c=[choice(pool)]
					del pool[pool.index(c[0])]
					c.append(atoms[c[0]].symbol)
					chosen.append(c)
					if c[1] not in unique_syms:
						unique_syms.append(c[1])
				
				if len(unique_syms) > 1:
					ready = True
			
			#print("hello1") 

			# Going through the current structure & working out where the viable vacancy sites are

			cell=current_structure['atoms'].cell.cellpar()
			steps=[grid_spacing/cell[0],grid_spacing/cell[1],grid_spacing/cell[2]]
			# initialize our grid
			grid=[] 
			
			#work out how many poitns we need along each axis
			#print(cell)
			npoints=[math.floor(cell[0]/grid_spacing),math.floor(cell[1]/grid_spacing),math.floor(cell[2]/grid_spacing)]
			#populate the grid
			
			for x in range(npoints[0]):
				for y in range(npoints[1]):
					for z in range(npoints[2]):
						grid.append([(x*steps[0])*cell[0],(y*steps[1])*cell[1],(z*steps[2])*cell[2]])
			
			#as a test, make a dummy cell with the grid in:
			full_grid=Atoms("",cell=current_structure['atoms'].cell,pbc=True)
			
			viable=[]
			
			for i in grid:
				at=Atom("X",position=i)
				full_grid+=at
				
				#now work out if the position is viable
				temp_atoms=current_structure['atoms'].copy()
				temp_atoms+=at
				distances=temp_atoms.get_distances(-1,list(range(len(temp_atoms)-1)),mic=True)
				if min(distances) > exclusion:
					viable.append(i)
			
			#print("viable vacancies: ",viable)
			#print("atoms chosen to swap: ",chosen)
			#print("hello2")
			#for each site to swap, go through and lable whether or not we're going to shuffle it, or try and move it to a vacancy
			# need to make sure that the number of sites > 2 or 0 and not the same element
			ready = False
			
			while ready == False:
				nsite=0
				nvac=0
				sites=[]
				vacs=[]
				for i in range(len(chosen)):
					c=random.random()
					# go with a 50:50 choice between the two
					if c < 0.5:
						sites.append(chosen[i])
						nsite+=1
						
					if c >= 0.5:
						vacs.append(chosen[i])
						nvac+=1
						
				if nsite == 0:
					ready = True
				if nsite >= 2:
					ready = True
			
			
			#print("hello3")
			
			#Go through and move atoms to vacancies:
			# note, on a one by one basis, as they get moved will have to update / check that each site is viable
			fail = False
			for i in vacs:
				done=False
				trials=0 # keep track of how many attempts we've made to move an atom, if we get stuck, switch to moving to a site.
				all_vacs=viable.copy() # take a copy of the vacancy site
				while done == False:
					if len(all_vacs)==0:
						fail=True
						break
						
					vac=choice(all_vacs)
					temp_atoms=atoms.copy()
					temp_atoms[i[0]].position=vac
					distances=temp_atoms.get_distances(i[0],list(range(len(temp_atoms))),mic=True)
					# if the move is within the exclusion radius, update the atoms object
					if min(distances) > exclusion:
						atoms=temp_atoms.copy()
						done=True
						
					else:
						#once we've tried and failed to use a vacancy, remove it from the pool
						del all_vacs[all_vacs.index(vac)]
						trials+=1
						
					if trials >= len(viable): # seems reasonable to attmpt each of the listed vacencies before giving up!
						sites.append(i)
						done=True
						
				if fail == True:
					break
					
			if fail == True:
				move= choice(moves_to_choose)
				continue			
					
				#first choose a vacency
			#print("hello1")

			
			#Go through and move the sites around:
			#build a new list with the symbols shuffled
			shuffled=[]
			indicies=[]
			symbols=[]
			
			for i in sites:
				indicies.append(i[0])
				symbols.append(i[1])
			
			o_symbols=symbols.copy()
			
			while o_symbols == symbols:
				shuffle(symbols)
			
			for i in range(len(indicies)):
				c=[indicies[i],symbols[i]]
				shuffled.append(c)
			
			# now go through and reasign the chemical symbols in the structure
			
			for i in shuffled:
				atoms[i[0]].symbol=i[1]
			
			
			# before we pass the structure back, check that we have a set of atoms which is consistent with the input composition
			trial_fus=len(atoms)/sum(list(composition.values()))
			counts=[]
			correct=False
			if float(trial_fus).is_integer():
				#now need to check each species
				symbols=atoms.get_chemical_symbols()
				for x in list(composition.keys()):
					num1=symbols.count(x)
					counts.append(num1/composition[x])
				if all(x == trial_fus for x in counts):
					correct = True
			
			if correct == True:
				comp_atoms=atoms.copy()
				complete=True

		##########################################################################
				
		##########################################################################

		# swapping the position of all atoms - including vacancies
		if move == 6:	
			# have started by first copying the code from move = 3: swapping some atoms,
			# then need to locate the vacancies within a structure & then decide how many
			# vacancy sites it wants to use when swapping atoms.
			atoms = current_structure['atoms'].copy()
			
			#check to see that there are enough atoms & species to swap in the structure
			
			syms=[]
			for i in atoms:
				if not i.symbol in syms:
					syms.append(i.symbol)
			
			# if we're dealing with an element, redraw a move and continue
			if len(syms) == 1:
				move = choice(moves_to_choose)
				continue			
			
			# we obviously need more than 2 atoms in the structure to make this move
			if len(atoms) <= 2:
				move = choice(moves_to_choose)
				continue			
			
			# then if we have enough atoms / species, get swapping!
			
			#flag to say we're ready to change the atoms objecct
			ready = False
			while ready == False:
				to_swap=len(atoms)				
				pool=list(range(len(atoms)))
				chosen=[[choice(pool)]]
				del pool[chosen[0][0]]
				unique_syms=[]
				#print(chosen[0])
				#print(atoms[chosen[0]])
				chosen[0].append(atoms[chosen[0]].symbols[0])
				for i in range(to_swap-1):
					c=[choice(pool)]
					del pool[pool.index(c[0])]
					c.append(atoms[c[0]].symbol)
					chosen.append(c)
					if c[1] not in unique_syms:
						unique_syms.append(c[1])
				
				if len(unique_syms) > 1:
					ready = True


			# Going through the current structure & working out where the viable vacancy sites are

			cell=current_structure['atoms'].cell.cellpar()
			steps=[grid_spacing/cell[0],grid_spacing/cell[1],grid_spacing/cell[2]]
			# initialize our grid
			grid=[] 
			
			#work out how many poitns we need along each axis
			#print(cell)
			npoints=[math.floor(cell[0]/grid_spacing),math.floor(cell[1]/grid_spacing),math.floor(cell[2]/grid_spacing)]
			#populate the grid
			
			for x in range(npoints[0]):
				for y in range(npoints[1]):
					for z in range(npoints[2]):
						grid.append([(x*steps[0])*cell[0],(y*steps[1])*cell[1],(z*steps[2])*cell[2]])
			
			#as a test, make a dummy cell with the grid in:
			full_grid=Atoms("",cell=current_structure['atoms'].cell,pbc=True)
			
			viable=[]
			
			for i in grid:
				at=Atom("X",position=i)
				full_grid+=at
				
				#now work out if the position is viable
				temp_atoms=current_structure['atoms'].copy()
				temp_atoms+=at
				distances=temp_atoms.get_distances(-1,list(range(len(temp_atoms)-1)),mic=True)
				if min(distances) > exclusion:
					viable.append(i)
			
			#print("viable vacancies: ",viable)
			#print("atoms chosen to swap: ",chosen)
			
			#for each site to swap, go through and lable whether or not we're going to shuffle it, or try and move it to a vacancy
			# need to make sure that the number of sites > 2 or 0 and not the same element
			ready = False
				
				
			while ready == False:
				nsite=0
				nvac=0
				sites=[]
				vacs=[]
				for i in range(len(chosen)):
					c=random.random()
					# go with a 50:50 choice between the two
					if c < 0.5:
						sites.append(chosen[i])
						nsite+=1
						
					if c >= 0.5:
						vacs.append(chosen[i])
						nvac+=1
						
				if nsite == 0:
					ready = True
				if nsite >= 2:
					ready = True
			
			
			#Go through and move atoms to vacancies:
			# note, on a one by one basis, as they get moved will have to update / check that each site is viable
			
			for i in vacs:
				done=False
				trials=0 # keep track of how many attempts we've made to move an atom, if we get stuck, switch to moving to a site.
				all_vacs=viable.copy() # take a copy of the vacancy site
				while done == False:
					try:		
						vac=choice(all_vacs)
					except:
						break
						
					temp_atoms=atoms.copy()
					temp_atoms[i[0]].position=vac
					distances=temp_atoms.get_distances(i[0],list(range(len(temp_atoms))),mic=True)
					# if the move is within the exclusion radius, update the atoms object
					if min(distances) > exclusion:
						atoms=temp_atoms.copy()
						done=True
						
					else:
						#once we've tried and failed to use a vacancy, remove it from the pool
						del all_vacs[all_vacs.index(vac)]
						trials+=1
						
					if trials >= len(viable): # seems reasonable to attmpt each of the listed vacencies before giving up!
						sites.append(i)
						done=True
				#first choose a vacency
			
			
			#Go through and move the sites around:
			#build a new list with the symbols shuffled
			shuffled=[]
			indicies=[]
			symbols=[]
			
			for i in sites:
				indicies.append(i[0])
				symbols.append(i[1])
			
			o_symbols=symbols.copy()
			#print("start symbols: ",o_symbols)
			while o_symbols == symbols:
				shuffle(symbols)
			
			for i in range(len(indicies)):
				c=[indicies[i],symbols[i]]
				shuffled.append(c)
			
			# now go through and reasign the chemical symbols in the structure
			
			for i in shuffled:
				atoms[i[0]].symbol=i[1]
			
			# before we pass the structure back, check that we have a set of atoms which is consistent with the input composition
			trial_fus=len(atoms)/sum(list(composition.values()))
			counts=[]
			correct=False
			if float(trial_fus).is_integer():
				#now need to check each species
				symbols=atoms.get_chemical_symbols()
				for x in list(composition.keys()):
					num1=symbols.count(x)
					counts.append(num1/composition[x])
				if all(x == trial_fus for x in counts):
					correct = True

			if correct == True:
				comp_atoms=atoms.copy()
				complete=True

		##########################################################################
				
		###########################################################################
		#if move == 7: # 7. swap the positions of two sub-modules
		#	#### ORIGNAL VERSION ####
		#	
		#	start_point=current_structure.copy()
		#	#get a copy of the input modules
		#	start_modules=start_point['modules']
		#	#chose the two sub modules that we want to swap
		#	c1=choice(list(range(len(start_modules))))
		#	c2=c1
		#	while c2 == c1:
		#		c2 = choice(list(range(len(start_modules))))
		#	
		#	#swap the position of the two modules
		#	
		#	#at the moment, the below doesn't work as it doesn't update the 3D position of the sub module. What I need to do, is work out the position in space of the origin of the two sub modules and shift them accordinly?
		#	
		#	swapped_modules=start_modules.copy()
		#	#
		#	swapped_modules[c1]=start_modules[c2]
		#	swapped_modules[c2]=start_modules[c1]
		#	
		#	#print(swapped_modules[c1])
		#	#print(swapped_modules[c1].get_positions())
		#	#print(swapped_modules[c1].cell.cellpar())
		#	
		#	
		#	#now need to go and rebuild the new atoms object
		#	
		#	t={'modules':swapped_modules,
		#		'sub module cell':start_point['sub module cell'],
		#		'shape in submods':start_point['shape in submods'],
		#		'ap':start_point['ap']}
		#		
		#	#put the atoms obeject together and error check it	
		#	
		#	ap=t['ap']
		#	target_number_atoms=len(start_point['atoms'])
		#	#check_bonds=False
		#	
		#	atoms=assemble_structure2(t)
		#	accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=target_number_atoms)
		#	trial+=1
		#	if accept == 1:
		#		complete=True
		#		
		#	#give up with this move if we can't do it after 5000 attempts
		#		
		#	if accept ==0:
		#		if trial > 5000:
		#			#trim moves to remove 7 then choose a different move:
		#			t_moves=[]
		#			for i in moves_to_choose:
		#				if i != 7:
		#					t_moves.append(i)
		#			
		#			move = choice(t_moves)
		#			trial = 0
		#			continue
	   #
		#	#write("before.cif",start_point['atoms'])
		#	#write("after.cif",atoms)
	   #
		#		
		###########################################################################

		##########################################################################
		if move == 7: # 7. swap the positions of two sub-modules, version taking into account the difference in where the atom closest to the origin is
			start_point=current_structure.copy()
			#write("temp.cif",start_point['atoms'])
			
			#get a copy of the input modules
			start_modules=start_point['modules'].copy()
			
			if len(start_modules) == 1:
				move=choice(moves_to_choose)
				continue
			
			#chose the two sub modules that we want to swap
			c1=choice(list(range(len(start_modules))))
			c2=c1
			while c2 == c1:
				c2 = choice(list(range(len(start_modules))))
			
			#swap the position of the two modules
			
			mod1=start_modules[c1].copy()
			mod2=start_modules[c2].copy()
			
			need_to_move=False
			
			if len(mod1) > 0:
				if len(mod2) > 0:
					need_to_move=True
					
			if need_to_move==True:
				#for each module, find the closest atom to the origin of the sub-module
				dum=Atom("X",position=[0,0,0])
				mod1+=dum
				mod2+=dum
		   	
				d1=list(mod1.get_distances(-1,list(range(len(mod1)-1)))) 
				m1=d1.index(min(d1))
         	
				d2=list(mod2.get_distances(-1,list(range(len(mod2)-1)))) 
				m2=d2.index(min(d2))
				
				#now need to work out the difference in fractional co-ordinates between the nearest atoms to the origin
				
				#first work out how fara we need to move the atoms in sub module 1
				diff1=mod1[m1].scaled_position - mod2[m2].scaled_position
				s_positions1=list(mod1.get_scaled_positions())
				del s_positions1[-1] # remove the last entry, as we don't care about the dummy atom
				#print(s_positions1)
				for w in range(len(s_positions1)):
					s_positions1[w]=s_positions1[w] - diff1
				
				# make the shifted sub module
				mod1s=start_modules[c1].copy()
				mod1s.set_scaled_positions(s_positions1)
				#view([mod1,mod1s])
				
							
				#now same again for sub module 2
				diff2=mod2[m2].scaled_position - mod1[m1].scaled_position
				s_positions2=list(mod2.get_scaled_positions())
				del s_positions2[-1]
				for w in range(len(s_positions2)):
					s_positions2[w]=s_positions2[w]-diff2
				
				mod2s=start_modules[c2].copy()
				mod2s.set_scaled_positions(s_positions2)
				#view([mod2,mod2s])
				
				#sys.exit()
				##at the moment, the below doesn't work as it doesn't update the 3D position of the sub module. What I need to do, is work out the position in space of the origin of the two sub modules and shift them accordinly?
				
				swapped_modules=start_modules.copy()
				#
				swapped_modules[c1]=mod2s.copy()
				swapped_modules[c2]=mod1s.copy()
			
			if need_to_move == False:
				swapped_modules=start_modules.copy()
				swapped_modules[c1]=start_modules[c2]
				swapped_modules[c2]=start_modules[c1]
			
			#print(swapped_modules[c1])
			#print(swapped_modules[c1].get_positions())
			#print(swapped_modules[c1].cell.cellpar())
			
			
			#now need to go and rebuild the new atoms object
			
			t={'modules':swapped_modules,
				'sub module cell':start_point['sub module cell'],
				'shape in submods':start_point['shape in submods'],
				'ap':start_point['ap']}
				
			#put the atoms obeject together and error check it	
			
			ap=t['ap']
			target_number_atoms=len(start_point['atoms'])
			#check_bonds=False
			
			atoms=assemble_structure2(t)
			accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=target_number_atoms)
			trial+=1
			if accept == 1:
				# before we pass the structure back, check that we have a set of atoms which is consistent with the input composition
				trial_fus=len(atoms)/sum(list(composition.values()))
				counts=[]
				correct=False
				if float(trial_fus).is_integer():
					#now need to check each species
					symbols=atoms.get_chemical_symbols()
					for x in list(composition.keys()):
						num1=symbols.count(x)
						counts.append(num1/composition[x])
					if all(x == trial_fus for x in counts):
						correct = True
				
				if correct == True:
					comp_atoms=atoms.copy()
					complete=True
				
			#give up with this move if we can't do it after 5000 attempts
			#write("before.cif",start_point['atoms'])
			#write("after.cif",atoms)
			#sys.exit()	
			if accept ==0:
				if trial > 100:
					#trim moves to remove 7 then choose a different move:
					t_moves=[]
					for i in moves_to_choose:
						if i != 7:
							t_moves.append(i)
					
					move = choice(t_moves)
					trial = 0
					continue
	
	
				
		##########################################################################


		###########################################################################
		if move == 8: 	#8. find two full slices/modules in the structure and switch their positions
			#### ORIGINAL VERSION ####
			start_point = current_structure.copy()
			shape=start_point['shape in submods'].copy()
			layers=shape[2]
			# If we have less than three layeers the move it pointless, move on!
			if layers < 3:
				t_moves=[]
				for i in moves_to_choose:
					if i != 8:
						t_moves.append(i)
				
				move = choice(t_moves)
				trial = 0
				continue
			
			# then group the sub modules into their full modules, so we can swap them around
			
			nsub=int(len(start_point['modules']))
			sub_mods=start_point['modules'].copy()
			
			#then work out the layers
			modules={}
			
			# work out how many sub modules we have in each layer
			size=shape[0]*shape[1]
			# how many modules we expect to have
			
			#now group the sub modules
			mods_used=0
			for i in range(layers):
				group=[]
				for j in range(size):
					group.append(sub_mods[0])
					del sub_mods[0]
				modules[i+1]=group
			
			#choose the layers to swap
			c1=choice(list(modules.keys()))
			c2=choice(list(modules.keys()))
			while c2 == c1:
				c2=choice(list(modules.keys()))
			
			#swap the keys around to swap the module order
			old_keys=list(modules.keys())
			new_keys=[]
			for i in range(len(old_keys)):
				if old_keys[i] == c1:
					new_keys.append(c2)
					continue
					
				if old_keys[i] == c2:
					new_keys.append(c1)
					continue
					
				else:
					new_keys.append(old_keys[i])
					continue
			
			#At the moment, this doesn't work as it doesn't update the 3D position of the atoms within the layer.
			#now go through and rebuild the 'modules' object for the altered structure
			sub_mods_2=[]
			for i in new_keys:
				mod=modules[i]
				for j in mod:
					sub_mods_2.append(j)
			
			#now build the temporary structure object to pass to assemble structure
			t={'modules':sub_mods_2.copy(),
				'sub module cell':start_point['sub module cell'].copy(),
				'shape in submods':start_point['shape in submods'].copy(),
				'ap':start_point['ap']}				
			ap=t['ap']
			target_number_atoms=len(current_structure['atoms'])
			#check_bonds=False
			try:
				atoms=assemble_structure2(t)
				accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=target_number_atoms)
			except:
				accept=0
			
			
			trial+=1
			
			if accept == 1: # before we pass the structure back, check that we have a set of atoms which is consistent with the input composition
				trial_fus=len(atoms)/sum(list(composition.values()))
				counts=[]
				correct=False
				if float(trial_fus).is_integer():
					#now need to check each species
					symbols=atoms.get_chemical_symbols()
					for x in list(composition.keys()):
						num1=symbols.count(x)
						counts.append(num1/composition[x])
					if all(x == trial_fus for x in counts):
						correct = True
				
				if correct == True:
					comp_atoms=atoms.copy()
					complete=True
						
			#give up with this move if we can't do it after 5000 attempts
				
			if accept ==0:
				if trial > 500:
					#trim moves to remove 8 then choose a different move:
					t_moves=[]
					for i in moves_to_choose:
						if i != 8:
							t_moves.append(i)
					
					move = choice(t_moves)
					trial = 0
					continue
							
		###########################################################################

		##########################################################################
		if move == 100: 	#8. find two full slices/modules in the structure and switch their positions version taking into account the difference in where the atom closest to the origin is
			start_point = current_structure.copy()
			shape=start_point['shape in submods']
			layers=shape[2]
			# If we have less than three layeers the move it pointless, move on!
			if layers < 3:
				t_moves=[]
				for i in moves_to_choose:
					if i != 8:
						t_moves.append(i)
				
				move = choice(t_moves)
				trial = 0
				continue
			
			# then group the sub modules into their full modules, so we can swap them around
			
			nsub=int(len(start_point['modules']))
			sub_mods=start_point['modules'].copy()
			
			#then work out the layers
			modules={}
			
			# work out how many sub modules we have in each layer
			size=shape[0]*shape[1]
			# how many modules we expect to have
			
			#now group the sub modules
			mods_used=0
			for i in range(layers):
				group=[]
				for j in range(size):
					group.append(sub_mods[0])
					del sub_mods[0]
				modules[i+1]=group
			
			#choose the layers to swap
			c1=choice(list(modules.keys()))
			c2=choice(list(modules.keys()))
			while c2 == c1:
				c2=choice(list(modules.keys()))
			
			#swap the keys around to swap the module order
			old_keys=list(modules.keys())
			new_keys=[]
			for i in range(len(old_keys)):
				if old_keys[i] == c1:
					new_keys.append(c2)
					continue
					
				if old_keys[i] == c2:
					new_keys.append(c1)
					continue
					
				else:
					new_keys.append(old_keys[i])
					continue		
			
			
			#print(old_keys)
			#print(new_keys)
			
			#need to build the two modules which are going to be swapped
			#build module 1
			layer_shape=start_point['shape in submods'].copy()
			layer_shape[2]=1
			
			l1={'modules':modules[c1],
				'sub module cell':start_point['sub module cell'],
				'shape in submods':layer_shape,
				'ap':start_point['ap']}				
			
			layer1=assemble_structure2(l1)

			#view(layer1)
			
			#build layer 2
			l2={'modules':modules[c2],
				'sub module cell':start_point['sub module cell'],
				'shape in submods':layer_shape,
				'ap':start_point['ap']}				
			
			layer2=assemble_structure2(l2)
			
			#view(layer2)
			need_to_move=False
			
			if len(layer1) > 0:
				if len(layer2) > 0:
					need_to_move = True
					
			if need_to_move == True:
				#now as for move 7, in each layer work out the displacement we want to enforce
				#for module 1:
				dum=Atom("X",position=[0,0,0])
				mod1=layer1.copy()
				mod2=layer2.copy()
				
				mod1+=dum
				mod2+=dum
				
				d1=list(mod1.get_distances(-1,list(range(len(mod1)-1)))) 
				m1=d1.index(min(d1))
         	
				d2=list(mod2.get_distances(-1,list(range(len(mod2)-1)))) 
				m2=d2.index(min(d2))
				
				#now need to work out the difference in fractional co-ordinates between the nearest atoms to the origin
				
				#first work out how fara we need to move the atoms in sub module 1
				diff1=mod1[m1].scaled_position - mod2[m2].scaled_position
				diff2=mod2[m2].scaled_position - mod1[m1].scaled_position
				
				#shift sub modules in layer 1
				for w in range(len(modules[c1])):
					mod=modules[c1][w].copy()
					s_pos=list(mod.get_scaled_positions())
					for x in range(len(s_pos)):
						s_pos[x]=s_pos[x] - diff1
					
					if len(mod) > 0:
						modules[c1][w].set_scaled_positions([s_pos])
				
				#again for layer 2
				for w in range(len(modules[c2])):
					mod=modules[c2][w].copy()
					s_pos=list(mod.get_scaled_positions())
					for x in range(len(s_pos)):
						s_pos[x]=s_pos[x] - diff2
					
					if len(mod) > 0:
						modules[c2][w].set_scaled_positions([s_pos])
								
								
			#now go through and rebuild the 'modules' object for the altered structure
			sub_mods_2=[]
			for i in new_keys:
				mod=modules[i]
				for j in mod:
					sub_mods_2.append(j)
			
			#now build the temporary structure object to pass to assemble structure
			t={'modules':sub_mods_2,
				'sub module cell':start_point['sub module cell'],
				'shape in submods':start_point['shape in submods'],
				'ap':start_point['ap']}				
			ap=t['ap']
			target_number_atoms=len(start_point['atoms'])
			#check_bonds=False
			try:
				atoms=assemble_structure2(t)
				accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=target_number_atoms)
			except:
				accept=0
			
			#write("before.cif",start_point['atoms'])
			#write("after.cif",atoms)
			#
			#sys.exit()
			
			trial+=1
			if accept == 1:
				# before we pass the structure back, check that we have a set of atoms which is consistent with the input composition
				trial_fus=len(atoms)/sum(list(composition.values()))
				counts=[]
				correct=False
				if float(trial_fus).is_integer():
					#now need to check each species
					symbols=atoms.get_chemical_symbols()
					for x in list(composition.keys()):
						num1=symbols.count(x)
						counts.append(num1/composition[x])
					if all(x == trial_fus for x in counts):
						correct = True
				
				if correct == True:
					comp_atoms=atoms.copy()
					complete=True
				
			#give up with this move if we can't do it after 5000 attempts
				
			if accept ==0:
				if trial > 500:
					#trim moves to remove 8 then choose a different move:
					t_moves=[]
					for i in moves_to_choose:
						if i != 8:
							t_moves.append(i)
					
					move = choice(t_moves)
					trial = 0
					continue
							
		##########################################################################


		##########################################################################

		if move == 9: #generate a new set of instructions for the current modules set; equiv here may be trying to work out different possible shapes, perhaps using the solutions function from fuse107
			start_point = current_structure.copy()
			start_shape=start_point['shape in submods']
			#print(start_shape)
			nsub=len(start_point['modules'])
			
			#create a list of viable unit cells for our module set
			viable=[]
			#print(cubic_solutions)
			
			for i in list(cubic_solutions.keys()):
				if i == nsub:
					for j in cubic_solutions[i]:
						temp_cell=[j,j,2*j,90,90,90]
						viable.append([j,j,2*j])
			
			for i in list(tetragonal_solutions.keys()):
				if i == nsub:
					for j in tetragonal_solutions[i]:
						viable.append([j[0],j[0],j[1],90,90,90])			

			for i in list(hexagonal_solutions.keys()):
				if i == nsub:
					for j in hexagonal_solutions[i]:
						viable.append([j[0],j[0],j[1],90,90,120])			

			for i in list(orthorhombic_solutions.keys()):
				if i == nsub:
					for j in orthorhombic_solutions[i]:
						viable.append([j[0],j[1],j[2],90,90,90])			

			for i in list(monoclinic_solutions.keys()):
				if i == nsub:
					for j in monoclinic_solutions[i]:
						#need to choose cell angles:
						viable.append([j[0],j[1],j[2],choice([60,75,90,105,120,150]),choice([60,75,90,105,120,150]),choice([60,75,90,105,120,150])]) 			

			try:
				new_cell=choice(viable)
			except:
					t_moves=[]
					for i in moves_to_choose:
						if i != 9:
							t_moves.append(i)
					
					move = choice(t_moves)
					trial = 0
					continue

			#print(new_cell)
			
			#set up the new structure object:
			t={'modules':start_point['modules'].copy(),
				'sub module cell':start_point['sub module cell'].copy(),
				'shape in submods':new_cell,
				'ap':start_point['ap']}				
			ap=t['ap']
			target_number_atoms=len(start_point['atoms'])
			
			try:
				atoms=assemble_structure2(t)
				accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=target_number_atoms)
				
			except:
				accept = 0
				continue
		
			trial+=1
			if accept == 1:
				
				# before we pass the structure back, check that we have a set of atoms which is consistent with the input composition
				trial_fus=len(atoms)/sum(list(composition.values()))
				counts=[]
				correct=False
				if float(trial_fus).is_integer():
					#now need to check each species
					symbols=atoms.get_chemical_symbols()
					for x in list(composition.keys()):
						num1=symbols.count(x)
						counts.append(num1/composition[x])
					if all(x == trial_fus for x in counts):
						correct = True
				
				if correct == True:
					comp_atoms=atoms.copy()
					complete=True
				
			#give up with this move if we can't do it after 5000 attempts
				
			if accept ==0:
				if trial > 500:
					#trim moves to remove 9 then choose a different move:
					t_moves=[]
					for i in moves_to_choose:
						if i != 9:
							t_moves.append(i)
					
					move = choice(t_moves)
					trial = 0
					continue
		
		##########################################################################

		##########################################################################
		if move == 11: #	Attmpt to grow the structure along an axis. here there's also the oppertunity to apply symmetry operations to the grown part of the cell, e.g. a mirror plane or inversion centre
			start_point=current_structure.copy()
			#print("start atoms: ",len(start_point['atoms']) )
			#check to see if the number of atoms is less than half of the limit
			if len(start_point['atoms']) > 0.5*max_atoms:
				t_moves=[]
				for i in moves_to_choose:
					if i != 11:
						t_moves.append(i)
				
				move = choice(t_moves)
				continue
			
			temp_atoms=start_point['atoms'].copy()
			#choose direction to expand
			axis=choice([0,1,2])
			rep=[1,1,1]
			rep[axis]=2
			
			#choose, just repeat, extend + translate?, extend + invert?, mirror? # at somepoint want to try and work out rotate!
			typ=choice(["repeat","translate","invert","mirror"])
			
			#typ="mirror"
			
			#write("before.cif",temp_atoms)
			
			if typ=="repeat":
				temp_atoms=start_point['atoms'].repeat(rep)
			
			if typ=="translate":
				
				#translate the extension atoms by 1/2 the unit cell in the plane perpendicular to the repeat direction
				trans=[temp_atoms.cell.cellpar()[0],temp_atoms.cell.cellpar()[1],temp_atoms.cell.cellpar()[2]]
				for i in range(len(trans)):
					if i == axis:
						pass
					else:
						trans[i]=trans[i]*0.5
				
				to_add=temp_atoms.copy()
				cell=temp_atoms.cell.cellpar()
				temp_atoms.cell=[cell[0]*rep[0],cell[1]*rep[1],cell[2]*rep[2],cell[3],cell[4],cell[5]]
				for i in to_add:
					at=i
					at.position += trans
					temp_atoms+=at
				
			if typ=="invert":
				#invert the fractional co-ordinates in the extended part of the cell, with the inversion point in the middle of the face connecting the cells
				#first expand the cell and set up the translation object
				to_add=temp_atoms.copy()
					
				cell=temp_atoms.cell.cellpar()
					
				temp_atoms.cell=[cell[0]*rep[0],cell[1]*rep[1],cell[2]*rep[2],cell[3],cell[4],cell[5]]
				
				rep3=[rep[0]-1,rep[1]-1,rep[2]-1]#once we've inverted co-ordinates, how far we need to move the atom
				rep3=[rep3[0]*cell[0],rep3[1]*cell[1],rep3[2]*cell[2]]
				
				rep2=[1,1,1]
				for i in range(len(rep)):
					if rep[i]==2:
						rep2[i]=0
				
				centre=[0.5,0.5,0.5]
				for i in range(len(rep2)):
					centre[i] = centre[i]*rep2[i]
				
				#print("inversion centre: ",centre)
				
				for i in to_add:
					at=i
					s=at.scaled_position
					#print("start: ",s)
					i_move=[2*(centre[0]-s[0]),2*(centre[1]-s[1]),2*(centre[0]-s[2])]
					new_pos=[s[0]+i_move[0],s[1]+i_move[1],s[2]+i_move[2]]
					for i in range(len(new_pos)):
						if new_pos[i] < 0:
							new_pos[i]=1+new_pos[i]
					#print("new_pos: ",new_pos)
					
					start_abs=[s[0]*cell[0],s[1]*cell[1],s[2]*cell[2]]
					new_abs  =[new_pos[0]*cell[0],new_pos[1]*cell[1],new_pos[2]*cell[2]]
					
					for i in range(len(rep3)):
						new_abs[i]=new_abs[i]+rep3[i]
					
					at.position=new_abs
					temp_atoms+=at
					
					#sys.exit()
				
			if typ=="mirror": #mirror the co-ordinates from the face we're extending through
				#invert the fractional co-ordinates in the extended part of the cell, with the inversion point in the middle of the face connecting the cells
				#first expand the cell and set up the translation object
				to_add=temp_atoms.copy()
					
				cell=temp_atoms.cell.cellpar()
					
				temp_atoms.cell=[cell[0]*rep[0],cell[1]*rep[1],cell[2]*rep[2],cell[3],cell[4],cell[5]]
				
				#rep3=[rep[0]-1,rep[1]-1,rep[2]-1]#once we've inverted co-ordinates, how far we need to move the atom
				#rep3=[rep3[0]*cell[0],rep3[1]*cell[1],rep3[2]*cell[2]]
								
				#print("inversion centre: ",centre)
				
				for i in to_add:
					at=i
					s=at.scaled_position
					s1=at.position
					#print("direction: ",axis)
					#print("start :",at.position)
					
					d=[0,0,0]
					d[axis]=2*(1-s[axis])
					d_abs=[s1[0]+(d[0]*cell[0]),s1[1]+(d[1]*cell[1]),s1[2]+(d[2]*cell[2])]
					#print("new: ",d_abs)
					at.position=d_abs
					temp_atoms+=at
			
					#sys.exit()
					
			#view(temp_atoms)
			
			#write("after.cif",temp_atoms)
				
			#print(trans)
			#sys.exit()
			#print("end atoms: ",len(temp_atoms))
			
			ap=start_point['ap']
			target_number_atoms=len(temp_atoms)
			accept=error_check_structure(temp_atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=target_number_atoms)
			
			trial+=1
			if accept == 1:
				atoms=temp_atoms.copy()
				
				# before we pass the structure back, check that we have a set of atoms which is consistent with the input composition
				trial_fus=len(atoms)/sum(list(composition.values()))
				counts=[]
				correct=False
				if float(trial_fus).is_integer():
					#now need to check each species
					symbols=atoms.get_chemical_symbols()
					for x in list(composition.keys()):
						num1=symbols.count(x)
						counts.append(num1/composition[x])
					if all(x == trial_fus for x in counts):
						correct = True
					
					
				if correct == True:	
					comp_atoms=atoms.copy()
					complete=True
				
			#give up with this move if we can't do it after 5000 attempts
				
			if accept ==0:
				if trial > 50:
					#trim moves to remove 11 then choose a different move:
					t_moves=[]
					for i in moves_to_choose:
						if i != 11:
							t_moves.append(i)
					
					move = choice(t_moves)
					trial = 0
					continue
			
			#complete=True
		##########################################################################
		
		##########################################################################
		if move == 12: # triple the length of a structure, akin to move 11, again with the oppertunity to apply symmetry operations to the new part, e.g. translating it by 1/3. each time
			
			start_point=current_structure.copy()
			#blocks=2
			#check to see if the number of atoms is less than 1/3 of the limit
			if len(start_point['atoms']) > float(1/3.)*max_atoms:
				t_moves=[]
				for i in moves_to_choose:
					if i != 12:
						t_moves.append(i)
				
				move = choice(t_moves)
				continue			
			
			temp_atoms=start_point['atoms'].copy()
			#choose direction to expand
			axis=choice([0,1,2])
			rep=[1,1,1]
			rep[axis]=3
			#print(rep)
			#choose, just repeat, extend + translate?,# at somepoint want to try and work out rotate!
			typ=choice(["repeat","translate"])			
			
			#typ="repeat"
			
			if typ=="repeat":
				temp_atoms=start_point['atoms'].repeat(rep)
					
			if typ=="translate":
					
				#translate the extension atoms by 1/2 the unit cell in the plane perpendicular to the repeat direction
				
				trans=[temp_atoms.cell.cellpar()[0],temp_atoms.cell.cellpar()[1],temp_atoms.cell.cellpar()[2]]
				trans2=[temp_atoms.cell.cellpar()[0],temp_atoms.cell.cellpar()[1],temp_atoms.cell.cellpar()[2]]
				
				for i in range(len(trans)):
					if i == axis:
						pass
					else:
						trans[i]=trans[i]*1/3.
				
				to_add=start_point['atoms'].copy()
				cell=temp_atoms.cell.cellpar()
				temp_atoms.cell=[cell[0]*rep[0],cell[1]*rep[1],cell[2]*rep[2],cell[3],cell[4],cell[5]]
				#block1
				
				for i in to_add:
					at=i
					at.position += trans
					temp_atoms+=at				
		   	
				to_add=start_point['atoms'].copy()
				#block2
				for i in range(len(trans)):
					if i == axis:
						trans[i] = trans2[i]*2
					else:
						trans[i]=trans2[i]*2/3.
         	
				for i in to_add:
					at=i
					at.position += trans
					temp_atoms+=at				
				
				#view(start_point['atoms'])
				#view(temp_atoms)
				
			ap=start_point['ap']
			target_number_atoms=len(temp_atoms)
			accept=error_check_structure(
			temp_atoms,
			ideal_density,
			density_cutoff,
			check_bonds,
			btol,
			system_type,
			fu,
			composition,
			bondtable,
			ap,
			check_distances,
			dist_cutoff,
			target_number_atoms=target_number_atoms
			)
			
			#write("before.cif",start_point['atoms'])
			#write("after.cif",temp_atoms)
			#print("accept: ",accept)
			
			#sys.exit()
			
			trial+=1
			if accept == 1:
				atoms=temp_atoms.copy()
				
				# before we pass the structure back, check that we have a set of atoms which is consistent with the input composition
				trial_fus=len(atoms)/sum(list(composition.values()))
				counts=[]
				correct=False
				if float(trial_fus).is_integer():
					#now need to check each species
					symbols=atoms.get_chemical_symbols()
					for x in list(composition.keys()):
						num1=symbols.count(x)
						counts.append(num1/composition[x])
					if all(x == trial_fus for x in counts):
						correct = True
				
				if correct == True:
					comp_atoms=atoms.copy()
					complete=True
				
			#give up with this move if we can't do it after 5000 attempts
				
			if accept ==0:
				if trial > 50:
					#trim moves to remove 12 then choose a different move:
					t_moves=[]
					for i in moves_to_choose:
						if i != 12:
							t_moves.append(i)
					
					move = choice(t_moves)
					trial = 0
					continue
				
				#sys.exit()
				
		##########################################################################
		
		##########################################################################
		#	13. random new structure with upto the same or fewer fus as we have currently		
		if move == 13:
			start_point=current_structure.copy()
			if len(start_point['atoms']) < 4: 
				move = 14
				continue
			#first need the number of formula units in the current structure
			fus=int(len(start_point['atoms'])/atoms_per_fu)
			ap=start_point['ap']
			
			#atoms_per_fu=sum(composition.values())
			#
			max_atoms=int(len(start_point['atoms']))
			imax_atoms=int(len(start_point['atoms']))
			#imax_atoms=max_atoms
			trial = 0
			while trial <= 50:
				new_structure = get_new_structure(
				composition=composition,
				max_atoms=max_atoms,
				imax_atoms=max_atoms,
				max_ax=max_ax,
				density_cutoff = density_cutoff,
				check_bonds= check_bonds,
				btol= btol,
				check_distances= check_distances,
				system_type= system_type,
				dist_cutoff = dist_cutoff,
				vac_ratio = vac_ratio,
				atoms_per_fu= atoms_per_fu,
				imax_fus= fus,
				max_fus= fus,
				cubic_solutions= cubic_solutions,
				tetragonal_solutions= tetragonal_solutions,
				hexagonal_solutions= hexagonal_solutions,
				orthorhombic_solutions= orthorhombic_solutions,
				monoclinic_solutions= monoclinic_solutions,
				bondtable=bondtable,
				ideal_density=ideal_density,
				fu=fu,
				ap=ap
				)				
						
				atoms=new_structure['atoms']
				trial+=1
				
			# before we pass the structure back, check that we have a set of atoms which is consistent with the input composition
			trial_fus=len(atoms)/sum(list(composition.values()))
			counts=[]
			correct=False
			if float(trial_fus).is_integer():
				#now need to check each species
				symbols=atoms.get_chemical_symbols()
				for x in list(composition.keys()):
					num1=symbols.count(x)
					counts.append(num1/composition[x])
				if all(x == trial_fus for x in counts):
					correct = True

			if correct == True:
				comp_atoms=atoms.copy()
				complete=True
				
			if accept ==0:
				if trial > 50:
					#trim moves to remove 12 then choose a different move:
					t_moves=[]
					for i in moves_to_choose:
						if i != 12:
							t_moves.append(i)
					
					move = choice(t_moves)
					trial = 0
					continue
			
		##########################################################################
		
		##########################################################################
		#	14. random new structure, allowing to pull from any remaining structures which have been provided
		if move == 14:
			start_point=current_structure.copy()
			
			#print("using_prebuilt: ",using_prebuilt)
			#print("pre_built_structures: ",pre_built_structures.keys())
			# work out if we still have structures to pull from
			pre_built_to_pull=False
			
			if using_prebuilt == True:
				#create a list of the keys which have been used
				available=[]	
				for i in list(pre_built_structures.keys()):
					if not pre_built_structures[i]['used?']==True:
						available.append(i)
				
				if len(available) > 0:
					pre_built_to_pull=True
				
			#print("pre built we can use?: ",pre_built_to_pull)
			
			#if we have prebuilt available, use 50:50 of using one ot the other
			if pre_built_to_pull == True:
				new=choice(["random","prebuilt"])
			
			#all else fails, choose random
			else:
				new="random"
							
			if new == "prebuilt":
				c=choice(available)
				temp=pre_built_structures[c]
				atoms=temp['structure']
				pre_built_structures[c]['used?']=True
				comp_atoms=atoms.copy()
				complete=True
			
			if new == "random":
				fus=int(len(start_point['atoms'])/sum(composition.values()))
				ap=start_point['ap']
				
				atoms_per_fu=sum(composition.values())
				
				imax_atoms=max_atoms
				
				new_structure = get_new_structure(
				composition=composition,
				max_atoms=max_atoms,
				imax_atoms=max_atoms,
				max_ax=max_ax,
				density_cutoff = density_cutoff,
				check_bonds= check_bonds,
				btol= btol,
				check_distances= check_distances,
				system_type= system_type,
				dist_cutoff = dist_cutoff,
				vac_ratio = vac_ratio,
				atoms_per_fu= atoms_per_fu,
				imax_fus= max_fus,
				max_fus= max_fus,
				cubic_solutions= cubic_solutions,
				tetragonal_solutions= tetragonal_solutions,
				hexagonal_solutions= hexagonal_solutions,
				orthorhombic_solutions= orthorhombic_solutions,
				monoclinic_solutions= monoclinic_solutions,
				bondtable=bondtable,
				ideal_density=ideal_density,
				fu=fu,
				ap=ap)				
							
				atoms=new_structure['atoms']
				
				# before we pass the structure back, check that we have a set of atoms which is consistent with the input composition
				trial_fus=len(atoms)/sum(list(composition.values()))
				counts=[]
				correct=False
				if float(trial_fus).is_integer():
					#now need to check each species
					symbols=atoms.get_chemical_symbols()
					for x in list(composition.keys()):
						num1=symbols.count(x)
						counts.append(num1/composition[x])
					if all(x == trial_fus for x in counts):
						correct = True
				
				
				if correct == True:
					comp_atoms=atoms.copy()
					complete=True
		##########################################################################
		
		##########################################################################
				
		##########################################################################
		
		if move == 10:# move inspired by genetic algorithms, e.g. mutating / mating the structure, here may not have to be a full cut between structures, may be able to locat sub-modules or modules that we can swap in
			start_point=current_structure.copy()
			#print(start_point.keys())
			opt=choice([1])
			#option 1, pick a sub module, redraw the positions to original FUSE sub-modules
			#in principle, can make an option 2, where we do this for either an entire layer or the whole structure?
			mods=start_point['modules'].copy()
			
			viable=[]
			
			for i in range(len(mods)):
				if len(mods[i]) > 0:
					if len(mods[i]) < 4:
						viable.append(i)
			
			if len (viable) == 0:
				accept=0
				t_moves=[]
				for i in moves_to_choose:
					if i != 10:
						t_moves.append(i)
				
				move = choice(t_moves)
				trial = 0
				continue
				
				
			m1=choice(viable)
			
			m2=start_point['modules'][m1].copy()
		
			nums=list(m2.numbers)
			
			posns=[[0,0,0],[0.5,0,0],[0,0.5,0],[0.5,0.5,0]] # original grid points
			
			missing=0
			if len(nums) < 4:
				missing=4-len(nums)
			
			for i in range(missing):
				nums.append(0)
				
			#print(posns)
			shuffle(nums)
			#print(nums)
			sets={'number':[],'position':[]}
			
			for i in range(len(nums)):
				if nums[i] > 0:
					sets['number'].append(nums[i])
					sets['position'].append(posns[i])
						
			for i in range(len(m2)):
				m2[i].number=sets['number'][i]
			
			m2.set_scaled_positions(sets['position'])
			
			temp=start_point.copy()
			temp['modules'][m1]=m2

			t={'modules':temp['modules'],
				'sub module cell':start_point['sub module cell'],
				'shape in submods':start_point['shape in submods'],
				'ap':start_point['ap']}
				
			#put the atoms obeject together and error check it	
			
			ap=t['ap']
			target_number_atoms=len(start_point['atoms'])
			#check_bonds=False
			
			atoms=assemble_structure2(t)

			# before we pass the structure back, check that we have a set of atoms which is consistent with the input composition
			trial_fus=len(atoms)/sum(list(composition.values()))
			counts=[]
			correct=False
			if float(trial_fus).is_integer():
				#now need to check each species
				symbols=atoms.get_chemical_symbols()
				for x in list(composition.keys()):
					num1=symbols.count(x)
					counts.append(num1/composition[x])
				if all(x == trial_fus for x in counts):
					correct = True
  			
  			#write a quick test, to see if the structure makes sense!!
			write("test.cif",atoms)
			try:
			    test=read("test.cif")
			    os.remove("test.cif")
			except:
			    correct = False
			    try:
			        os.remove("test.cif")
			    except:
			        pass
  			
  			
  			
			if correct == True:
				comp_atoms=atoms.copy()
				complete=True
			
			trial +=1
			if trial > 500:
				#trim moves to remove 12 then choose a different move:
				t_moves=[]
				for i in moves_to_choose:
					if i != 10:
						t_moves.append(i)
				
				move = choice(t_moves)
				trial = 0
				continue
		
		
		##########################################################################
			
			#view(atoms)	
			
			#sys.exit()
		
	#for testing, write out the structure with completed move		
	write("swapped.cif",atoms)
	atoms=comp_atoms.copy()
	#ok, now we've completed the move we need to fill out the structure object:
	if use_spglib == True:
		import spglib
		#from ase import Atoms
		try:
			#print("fine till spglib")
			#write("spglib_fault.cif",atoms)
			lattice,positions,numbers=spglib.standardize_cell(atoms,symprec=1.e-5)
			temp2=Atoms(numbers=numbers,pbc=True)
			temp2.cell=lattice
			#temp2.numbers=numbers		
			temp2.set_scaled_positions(positions)
			atoms=temp2.copy()
			
		except:
			# if spglib cannont find a new cell, it returns a Nonetype error, so we can just pass & keep the current atoms object
			pass
	
	#print("fine after spglib")
	
	## before we pass the structure back, check that we have a set of atoms which is consistent with the input composition
	#trial_fus=len(atoms)/sum(list(composition.values()))
	#counts=[]
	#correct=False
	#if float(trial_fus).is_integer():
	#	#now need to check each species
	#	symbols=atoms.get_chemical_symbols()
	#	for x in list(composition.keys()):
	#		num1=symbols.count(x)
	#		counts.append(num1/composition[x])
	#	if all(x == trial_fus for x in counts):
	#		correct = True
	#sys.exit()
	
	#print("swapped atoms: ",len(atoms))

	write("temp.cif",atoms)
	structure=extract_module(["temp.cif"],bondtable)
	structure['optimised?']=False
	structure['energy']=0.0
	structure['converged']=False
	os.remove("temp.cif")
	# just check we've got everything
	
	#print("'modules'          ", structure['modules']          )
	#print("'sub module cell'  ", structure['sub module cell']  )
	#print("'shape in submods' ", structure['shape in submods'] )
	#print("'nmods'            ", structure['nmods']            )
	#print("'ap'               ", structure['ap']               )
	#print("'atoms'            ", structure['atoms']            )
	#print("'optimised?'       ", structure['optimised?']       )
	#print("'energy'           ", structure['energy']           )
	#print("'converged'        ", structure['converged']        )
	
	return structure,move,pre_built_structures