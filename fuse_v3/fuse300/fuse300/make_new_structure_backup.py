#fuse imports
from fuse300.possible_solutions import *
from fuse300.generate_random_structure import *
from fuse300.extract_modules import *
import fuse300.bond_table

#other imports
from random import choice
import math
import numpy
from ase.io import *
from ase.visualize import *
import pandas
import spglib
import time
from ase.visualize import *
from random import shuffle
import os,sys,glob,random
from ase import Atoms
#from fuse300.run_chgnet import *

def get_new_structure(composition='',max_atoms='',imax_atoms='',restart='',
	max_ax='',density_cutoff ='',check_bonds='',btol='',check_distances='',
	system_type='', dist_cutoff = '',vac_ratio = '',atoms_per_fu='',imax_fus='',
	max_fus='',cubic_solutions='',tetragonal_solutions='',hexagonal_solutions='',
	orthorhombic_solutions='',monoclinic_solutions='',bondtable='',
	ideal_density='',fu='',ap='',use_spglib='',
	old=False, # use this True / False flag to switch between old / new algorithmic generators
	# options for new generator:
	start_box=4.0, # starting box size
	ratt_dist=0.05, # how much to randomly rattle the box
	run_calc=False, # run optimisations with chgnet as the structure is built?
	run_final_calc=False, # run optimisation with chgnet when the structure is completed?
	
):
		
	if use_spglib == True:
		import spglib

	if old == True:
		trial=0
		gen = False
		while gen == False:
			
			##############################################################################
			# now call the generate random structure function
			#for each iteration of generating n structures:
			target_fu=choice(list(range(1,imax_fus+1)))
			target_atoms=target_fu*atoms_per_fu
			#print("hello")
			#print(target_fu)
    	
			atoms,string,instructions,accept=generate_random_structure(target_atoms,cubic_solutions,tetragonal_solutions, hexagonal_solutions, orthorhombic_solutions, monoclinic_solutions, atoms_per_fu, ideal_density, density_cutoff, check_bonds, btol, system_type, fu, composition, bondtable, ap, check_distances, dist_cutoff, vac_ratio=vac_ratio, max_fus=imax_fus,target_number_atoms=target_atoms)
    	
			##############################################################################
			# now call write the atoms object to a temporary file & then call the extract modules function
				#from ase import Atoms
			if use_spglib==True:
				try:
					temp=(atoms.cell[:],atoms.get_scaled_positions(),atoms.get_atomic_numbers())
					lattice,positions,numbers=spglib.standardize_cell(temp,symprec=1.e-5)
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
    	
			if trial >= 100:
				structure=None
				break
    	
			trial+=1
			
		
    	
		#print(structure['sub module cell'])
		#print(structure['shape in submods'])
    	
		return structure
		
	if old == False:
		form = composition.copy()
		nform=list(range(1,max_fus+1))
		#print (f"form: {form}")
		#print (f"nform: {nform}")
		#print (f"start box: {start_box}")
		#print (f"ratt_dist: {ratt_dist}")
		#print (f"short: {dist_cutoff}")
		#print (f"run calc: {run_calc}")
		#print (f"run final calc: {run_final_calc}")
		#sys.exit()
		#CHGNET rough chgnet options
		n_opts=1
		rel=''
		relaxer_opts={
		'fmax':[2.0,],
		'steps':[500],
		'verbose':[False]
		}
		opt_class=['FIRE']
		opt_device='cpu' #Device to use for chgnet optimisation, 'cpu' or 'cuda'
		mode='relax' # set the calculation mode for chgnet, other option is 'single'
		use_spglib=True
		short=1.5
		n=choice(nform)
		formula=form.copy()
		for i in list(formula.keys()):
			formula[i]=form[i]*n

		#if not os.path.isdir("steps"):
		#	os.mkdir("steps")
		
		lattice=choice([1,2,3,4,5,6,7]) # choose the lattice we want to use, # 1=triclinic, 2=monoclinic, 3=orthorhombic, 4=tetragonal, 5=trigonal, 6=hexagonal, 7=cubic
		#Let's start with triclinic
		lattice = 1
		#space group format: 'group_symbol':[number,{sites<mult:coords>},[list of fixed_positions]]
		if lattice==1: # triclinic
			cell=[1*start_box,random.uniform(0.1,0.99)*start_box,random.uniform(0.1,0.99)*start_box,choice(list(range(20,120)))+random.random(),choice(list(range(20,120)))+random.random(),choice(list(range(20,120)))+random.random()]
		   #groups={
		   #'P1':[1,{1:[['x','y','z']]},
			#[],
			#],
			#'P-1':[2,{2:[['x','y','z']],1:[[0.5,0.5,0.5],[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0],[0.5,0,0],[0,0.5,0],[0,0,0.5],[0,0,0]]},
			#[[0.5,0.5,0.5],[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0],[0.5,0,0],[0,0.5,0],[0,0,0.5],[0,0,0]]
			#],
			#}
         #
			#sg=choice(['P1'])
		
		place=[]
		
		for i in list(formula.keys()):
			for j in range(formula[i]):
				place.append(i)
				
		#print(place)
		
		shuffle(place)
		
		ini_box=Atoms(cell=[start_box,start_box,start_box],pbc=True)
		a1=Atoms(place[0],positions=[[0,0,0]])
		ini_box+=a1
		del place[0]
		ini_box[0].scaled_position=[[0.5,0.5,0.5]]
		
		#view(ini_box)
		
		#atoms,energy,converged=run_chgnet(ini_box)
		
		if run_calc==True:
			atoms,energy,converged=run_chgnet(ini_box,n_opts=n_opts,rel=rel,relaxer_opts=relaxer_opts,opt_class=opt_class,opt_device=opt_device,use_spglib=False,mode=mode)
		else:
			atoms=ini_box.copy()
		
		#write("ini.cif",atoms)

		step=1
		
		while len(place) > 0:
			placed = False
			#temp_atoms=atoms.copy()
			#temp_atoms.cell=temp_atoms.cell*1.25
			t=0
			
			while placed == False:
				temp_atoms=atoms.copy()
				
				#print(t)
				if t == 10:
					original_positions=temp_atoms.get_scaled_positions()
					temp_atoms.cell=temp_atoms.cell*1.25
					temp_atoms.set_scaled_positions(original_positions)
				if t == 20:
					original_positions=temp_atoms.get_scaled_positions()
					temp_atoms.cell=temp_atoms.cell*1.5
					temp_atoms.set_scaled_positions(original_positions)
    	
				if t == 30:
					original_positions=temp_atoms.get_scaled_positions()
					temp_atoms.cell=temp_atoms.cell*1.75
					temp_atoms.set_scaled_positions(original_positions)
    	
				if t == 40:
					original_positions=temp_atoms.get_scaled_positions()
					temp_atoms.cell=temp_atoms.cell*2
					temp_atoms.set_scaled_positions(original_positions)
    	
    	
				if use_spglib == True:
					init_length=len(atoms)
					try:
						temp=(temp_atoms.cell[:],temp_atoms.get_scaled_positions(),temp_atoms.get_atomic_numbers())
						lattice,positions,numbers=spglib.standardize_cell(temp,symprec=1.e-5)
						temp2=Atoms(numbers=numbers,pbc=True)
						temp2.cell=lattice
						temp2.set_scaled_positions(positions)
						if len(temp2) == init_length:
							temp_atoms=temp2.copy()
						#print("I'm using SPGLIB!")
		
					except:
						#print("I failed at using SPGLIB!!")
						pass
		
				
				try:
					new_at=Atoms(place[0],positions=[[0,0,0]])
				except:
					pass
					
				try:
					temp_atoms+=new_at
				except:
					pass
				
				try:
					temp_atoms[-1].scaled_position=[[random.random(),random.random(),random.random()]]
				except:
					pass
				
				try:
					dists=temp_atoms.get_all_distances(mic=True)
				except:
					pass
					
				try:
					if min(dists[dists !=0]) >= short:
						placed=True
					#if min(dists[dists !=0]) <= short:
					#	print("oh bother")
				except:
					pass
					
				t+=1
				if t == 50:
					#temp_atoms.cell=temp_atoms.cell*1.25
					break
					
			if placed == True:
				#atoms,energy,converged = run_chgnet(atoms,n_opts=n_opts,rel=rel,relaxer_opts=relaxer_opts,opt_class=opt_class,opt_device=opt_device,use_spglib=use_spglib,mode=mode)
				start_atoms=temp_atoms.copy()
				if run_calc == True:
					temp_atoms,energy,converged=run_chgnet(temp_atoms,n_opts=n_opts,rel=rel,relaxer_opts=relaxer_opts,opt_class=opt_class,opt_device=opt_device,use_spglib=False,mode=mode)
				dists=temp_atoms.get_all_distances(mic=True)
				if min(dists[dists !=0]) <= short:
					temp_atoms=start_atoms.copy()
    	
				#temp_atoms,energy,converged=run_chgnet(temp_atoms,n_opts=n_opts)
				atoms=temp_atoms.copy()
				del place[0]
				#write(f"steps/step_{step}.cif",atoms)
				step+=1
				#break
		
		if run_final_calc == True:

			#CHGNET more refined options!
			n_opts=3
			rel=''
			relaxer_opts={
			'fmax':[1.0,0.5,0.05],
			'steps':[250,750,750],
			'verbose':[False,False,False]
			}
			opt_class=['FIRE','BFGS','BFGS']
			opt_device='cpu' #Device to use for chgnet optimisation, 'cpu' or 'cuda'
			mode='relax' # set the calculation mode for chgnet, other option is 'single'
			atoms,energy,converged=run_chgnet(atoms,n_opts=n_opts,rel=rel,relaxer_opts=relaxer_opts,opt_class=opt_class,opt_device=opt_device,use_spglib=False,mode=mode)

		write("temp.cif",atoms)
		input_files=["temp.cif"]
		structure=extract_module(input_files,bondtable)
		return structure

		
		#sys.exit()
