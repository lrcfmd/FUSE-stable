################################################################################
#fuse imports
################################################################################
from fuse202.possible_solutions import *
from fuse202.generate_random_structure import *
from fuse202.extract_modules import *
import fuse202.bond_table
from fuse202.make_new_structure import *
from fuse202.run_gulp import *
from fuse202.run_vasp import *
from fuse202.run_qe import *
from fuse202.run_chgnet import *
from fuse202.run_multiple_calculators import run_calculators
from fuse202.make_basin_move import *
from fuse202.plot_results import *
from fuse202.assemble_spp import *

################################################################################
#other imports
################################################################################
from random import choice
import math
import numpy
from ase.io import *
from ase.visualize import *
import pandas
import os
import sys
import glob
from decimal import *
import shutil
import random
import time
import warnings
import datetime
import pickle
from decimal import *
from ase import Atoms

################################################################################

#warnings.filterwarnings("ignore") # currently use this as python raises RuntimeError 
# warnings when a physically unreasonable unit cell is generated, this is caught
# and structures rejected by the "converged" variable, so the warnings aren't really
# needed. Comment this out for debugging.

t1=datetime.datetime.now()
#desctription of each of the BH moves to print out
#	1. swap two atoms - DONE
#	2. swap an atom in to a vacent space - DONE

#	3. swapping the position of more than two atoms - no vacancies - DONE
#	4. swapping the position of more than two atoms - including vacancies - DONE
#	5. swap all of the atoms - no vacancies -DONE
#	6. swap all of the atoms - including vacancies -DONE

#  Think that for both 7 & 8, need to think about lattice translations in the modules in order to avoid atom crashes, at the moment, we seem to be getting alot of atom crashes. other option, when creating the modules, try to flatten them to get them back to where we were with fuse originally?
#	7. swap the positions of two sub-modules - DONE
#	8. find two full slices/modules in the structure and switch their positions - DONE

#	9. generate a new set of instructions for the current modules set; equiv here may be trying to work out different possible shapes, perhaps using the solutions function from fuse107 - DONE

#	10. move inspired by genetic algorithms, e.g. mutating / mating the structure, here may not have to be a full cut between structures, may be able to locat sub-modules or modules that we can swap in

#	11. attmpt to grow the structure along an axis. here there's also the oppertunity to apply symmetry operations to the grown part of the cell, e.g. a mirror plane or inversion centre - DONE
#  12. triple the length of a structure, akin to move 11, again with the oppertunity to apply symmetry operations to the new part, e.g. translating it by 1/3. each time - DONE

#	13. random new structure with upto the same number of fus as we have currently - DONE
#	14. random new structure, allowing the use of any remaining pre-built structures - DONE

# KNOWN ISSUES: when iterations = 0, the code should just read any restart files and / or generate the initial population and exit. at the moment, if the relaxtaions
# for the inintial population are not complete, it will continue to run through them instead! Also need to get it to flush the output after each structure during
# the initial population.
# Another new bug... when getting a structure from VASP, if for some reason vasp doesn't appear to run corrcetly, it's carrying through the energy from SPP instead?!?
# Another thing to do: keep an internal log of the total walltime used for a simulation!

move_des={
1:"Swap two atoms",
2:"Swap an atom in to a vacent space",
3:"Swapping the position of more than two atoms",
4:"Swapping the position of more than two atoms - including vacancies",
5:"Swap all of the atoms",
6:"Swap all of the atoms - including vacancies",
7:"Swap the positions of two sub-modules",
8:"Swap the positions of two modules",
9:"Generate a new unit cell for the current module set",
10:"Mutation of structure inspired by Genetic Algorithms",
11:"Double the structure",
12:"Triple the structure",
13:"Random new structure upto current size",
14:"Random new structure, including pulling from prebuilt pool"}
################################################################################
#define main fuse function
################################################################################
def run_fuse(
composition= '', # composition in dictionary format, set the keys to be the element types and the values to be the composition in integers, e.g. {'Ca':3,'Al':2,'Si':3,'O':12} 
max_atoms= '', # maximum number of atoms to use within the calculation, e.g. 160
imax_atoms= '', # maximum number of atoms to use when generating a random initial population, e.g. 80
restart=False, # restart previous calculation
max_ax=40, # maximum number of sub-modules along a given unit cell axis
density_cutoff = 0.4, # density cutoff for randomly generated structures, expressed as a fraction of the maximum packing density (which FUSE works out internally) 
check_bonds=True, # check bond numbers / distances as part of the structure error checking
btol=0.25, # fraction of incorrect bonds that FUSE will accept when generating / modifying structures
check_distances=True, # check for short interatomic distances when error checking a structure?
dist_cutoff = 1.0, # shortest permitted interatomic contact defined in angstroms
system_type="neutral", # the system type? note: only using "neutral" at the moment!
vac_ratio = 4, # when generating random structures, the maximum ratio of vacent sites : 1 atom
write_all_structures=True, # write out all of the geometry optimised structures which FUSE generates
swap_searches = False, # when set to FUSE, allow it to switch between search routines
search = 1, # search routine for FUSE to use, 1 = basin hopping, 2 basin hopping with reinforcement learning, 3 = genetic algorithm
ap_scale='', # scale factor which can be applied to the FUSE calculated lattice parameter for sub modules
read_exisiting_structures = False, # read in and slice previously generated structures for modules?
path_to_structures = '', # should be a directory containing the structures you want to read in
initial_gen = 20, # number of structures which should be in the initial population
iterations = 0, # number of structures to relax in this run of the code
ratt_dist=0.05, # perturb structures around a standard deviation when they have been generated in Angstroms
T=0.02, #temperature parameter to use for basin hopping, recommend 0.005 for DFT calculations as we tend to get much smaller energy differences
e_prec=1.e-5, #round computed energies to this value, this can avoid the BH search bouncing around becuase it's hitting the same structure multiple times but with very small changes in energy.q
fixed_ap='', #allow the user to specify the ap value in the input file, e.g. output from an ML model 
output_graph_at_end=True, # use matplot lib to create a plot of the csp run at the end?
revert_atoms=None, # whether or not to revert back to the starting structure, or continue with the structure "as relaxed" after each step, set to False to continue with structures
pull_random=False, # if set to True, when using pre-built structure pull them in a random order
pull_spp_rank=True, # if True, pull pre-built structures as ranked by spp potentials. BEWARE! you need to first run a script to rank all pre-generated structures & produce one or more csv files containing file names and corresponding energies
use_spglib=True, # if True, everytime an atoms object is generated, tidy up the structure with spglib before the modules are extracted 

#variables used to run ML structure generation:
#gn boss model for structure generation
generate_gn_boss_structures=False, # if set to true, when FUSE is firt launched, it will run gn-boss to generate the pool of referennce structures for this calculation.
gn_boss_command=r"C:\ProgramData\anaconda3\_conda.exe run -p C:\ProgramData\anaconda3\envs\ML_FUSE python .\get_cifs_for_FUSE.py", # For my machine, I've setup gn-boss in a seperate python environment, this is the command for that version of python
gn_search='tpe', # 'rand' random search, 'tpe' baysian opt, 'pso' particle swarm
gn_max_step=500, # number of generation attempts for gn-boss
gn_template_path=r'C:\Users\cc0u5\Documents\csp_installs\ML-FUSE\template', #path to template files for using gn-boss 
gn_zn_range=[1,4], # numbers of formula units to scan with GN-BOSS for generating structures
rank_gn_structures='single', #if None; do not rank structures, this should only be set if pull_random = True above, if 'opti' rank with SPPs AFTER geometry optimising the, if 'sing' rank based on single point calculations with SPPs. 
clear_previous_gn_structures=False, #if set to True, before starting the calcluation, remove any previous structures from reference structures & gn-boss generated_results.
generate_structures_only=False, #if set to true, run the gn structure generation then exit, useful to pre-generate structures, then srtart the calculation properly with generate_gn_boss_structures = False
ranking='mixed', # option for what to use to rank gn-boss output, can be set to either "gulp" or "chgnet"
r_calcs=['gulp','chgnet'],

#gulp options to rank spps:
r_kwds=['opti conj conp noelectro','opti conj conp noelectro','opti conp noelectro','sing conp noelectro'], # keywords for the gulp input, 
#multiple strings indicate running GULP in multiple stages
r_gulp_opts=[
['\ninclude ./lib.lib\ndump temp.res\nmaxcyc 250\nstepmax 0.001\ntime 5 minutes\ngmax 0.1\ngtol 0.1\nftol 0.1\nxtol 0.1'],
['\ninclude ./lib.lib\ndump temp.res\nmaxcyc 250\nstepmax 0.5\ntime 5 minutes\ngmax 0.1\ngtol 0.1\nftol 0.1\nxtol 0.1'],
['\ninclude ./lib.lib\ndump temp.res\nmaxcyc 1500\nlbfgs_order 5000\ntime 5 minutes\ngmax 0.1\ngtol 0.1\nftol 0.1\nxtol 0.1'],
['\ninclude ./lib.lib\ndump temp.res\nmaxcyc 1500\nlbfgs_order 5000\ntime 5 minutes\ngmax 0.1\ngtol 0.1\nftol 0.1\nxtol 0.1'],
],	# options for gulp, must include one line per set of inputs in "kwds"
r_lib='dummy.lib', # library file for interatomic potentials				



#variables which need defining for each of the search rountines
#1. basin hopping
search_gen_bh=1,
n_moves={1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1,12:1,13:1,14:1}, # moves used during the normal part of the search
r_moves={1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1,12:1,13:1,14:1}, # moves used when attempting to escape a local minima
melt_threshold=150, # number of structures to wait for a downhil move before melting the system
rmax=1000, # convergence criteria
grid_spacing=0.75, # grid spacing to use in Angstroms for basin hopping moves when trying to find gaps in a structure
exclusion=1.5, # when moving atoms around, set minimum distance a potential site needs to be from an exisiting atom in Angstroms
#2. additional parameters for using RL
#3. GA

#variables which need defining for the calculators
ctype='', # calculator type to use: "gulp", "vasp", "qe", "mixed" (qe = Quantum Espresso)
#GULP:
kwds='', # keywords needed for gulp input
gulp_opts='', # options needed for gulp
lib='', # library file needed for gulp
shel='', # any shell species needed for gulp
gulp_command='gulp < gulp.gin > gulp.got', # system command used to run gulp
gulp_timeout='', # timeout command for running GULP, note only works in Windows!
calcs='',# list to tell FUSE which calculators to use and in which order
assemble_spp_=False, # if set to True, collate the spp potential library for the calculation
spp_path=r"C:\Users\cc0u5\Documents\SPP_library\SPP", #path to spp potential libraries
#VASP:
vasp_opts='', # options for each Vasp calculation
kcut=30, # parameter for number of kpoints, higher number = more kpoints, just set to 1 if using KSPACING flag in VASP
#QE
qe_opts='', # options for Quantum Espresso

#CHGNET
n_opts=2,
rel=StructOptimizer(),
relaxer_opts={
'fmax':[0.1,0.05],
'steps':[250,750]
},
opt_class=['FIRE','BFGSLineSearch'],
opt_device='cpu', #Device to use for chgnet optimisation, 'cpu' or 'cuda'
):
	t0=time.time()
################################################################################

	print("################################################################################")
	print("#									       #")
	print("#		       Flexible Unit Structure Engine			       #")
	print("#				 (FUSE)					       #")
	print("#				 v2.02            			       #")
	print("#				                          		       #")
	print("################################################################################")
	print("\n\n")
	
	#############################################################################
	# start up bits / calculations
	#############################################################################
	
	# processing parts of the input file
	# number of structures we have performed geometry optimisation on during this run
	itr=0 
	# get atoms per formula unit
	atoms_per_fu=sum(composition.values())
	# number of formula units in the initial population
	imax_fus=int(imax_atoms/atoms_per_fu)
	# number of formula units in the general population
	max_fus=int(max_atoms/atoms_per_fu)
	# get possible unit cell sizes / shapes
	cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,monoclinic_solutions=possible_solutions(max_ax,restart)
	# load bond table data
	bondtable=numpy.load("bondtable.npz",allow_pickle=True)
	bondtable=bondtable['bond_table'].item()
	# set the original temperature value
	T0=T
	#set the environment variable for spp_path if needed
	#if spp_path != None:
	#	os.environ['SPP_PATH']=spp_path
	
	#if needed import spglib
	if use_spglib==True:
		#print("I'm using SPGLIB!!!")
		import spglib
	
	#Check to see if SPP needs assembling:
	if assemble_spp_==True:
		elements=list(composition.keys())
		assemble_spp(elements,spp_path=spp_path)
	
	if max_fus == 0:
		print ("ERROR: maximum number of Fus equals 0!, please increase maximum number of atoms in input!")
		sys.exit()
	
	if imax_fus == 0:
		imax_fus=1	
		
	#calculating ideal density###################################################
	
	# create a string containing the atomic numbers for 1 FU #
	keys=list(composition.keys())
	fu=[]
	for i in range(len(keys)):
	    for j in range(composition[keys[i]]):
	        fu.append(Atoms(keys[i]).get_atomic_numbers()[0])
	        
	mass=0
	volume=0
	
	# calculate the total volume / mass of 1 FU #	
	for i in range(len(fu)):
	    temp=Atoms(numbers=[fu[i]])
	    mass+=temp.get_masses()[0]
	    temp=temp.get_chemical_symbols()[0]
	    temp=bondtable[temp]
	    if len(list(temp.keys())) > 1:
	        rs=[]
	        keys=list(temp.keys())
	        for j in range(len(list(temp.keys()))):
	            rs.append(temp[keys[j]][-1])
	        volume+=((4/3)*math.pi*(min(rs)**3))
	    else:
	        volume+=((4/3)*math.pi*(temp[list(temp.keys())[0]][-1]**3))
	
	ideal_density=mass/volume
	
	#############################################################################
	
	### compute ap value ########################################################
	
	#### compute the ap to be used for the sub-modules ##########################
	# convert the fu back to symbols
	symbol_fu=[]
	for i in range(len(fu)):
	    temp=Atoms(numbers=[fu[i]])
	    symbol_fu.append(temp.get_chemical_symbols()[0])
	# read in shannon radi table
	#temp=numpy.load("bondtable.npz",allow_pickle=True)
	#bond_table=temp['bond_table'].item()
	ap=0
	if system_type=="neutral":
	    for i in range(len(symbol_fu)):
	        average = 0
	        temp=list(bondtable[symbol_fu[i]].keys())
	        for j in range(len(temp)):
	            average+=bondtable[symbol_fu[i]][temp[j]][-1]
	        average=average/len(temp)
	        ap+=average
	    
	ap=ap/len(symbol_fu)
	ap=4*ap
	
	if ap_scale != '':
		ap=ap*ap_scale
	
	if fixed_ap != '':
		ap=fixed_ap
	
	ap = float(Decimal(ap).quantize(Decimal('1e-4')))
	#############################################################################
	
	#### work out the normalised version of the input composition ############### 
	norm_comp={}
	norm_factor=round( min(list(composition.values())) , 2)
	for i in list(composition.keys()):
		norm_comp[i]=round( composition[i]/norm_factor,2 )
		
	#############################################################################
	#############################################################################

	## test call to obtaining a new randomly generated structure in the new data format:
	#structure = get_new_structure(composition=composition,
	#max_atoms=max_atoms,imax_atoms=imax_atoms,restart=restart,max_ax=max_ax,
	#density_cutoff = density_cutoff,check_bonds= check_bonds,btol= btol,
	#check_distances= check_distances,system_type= system_type,
	#dist_cutoff = dist_cutoff,vac_ratio = vac_ratio,atoms_per_fu= atoms_per_fu,
	#imax_fus= imax_fus,max_fus= max_fus,cubic_solutions= cubic_solutions,
	#tetragonal_solutions= tetragonal_solutions,
	#hexagonal_solutions= hexagonal_solutions,
	#orthorhombic_solutions= orthorhombic_solutions,
	#monoclinic_solutions= monoclinic_solutions,bondtable=bondtable,
	#ideal_density=ideal_density,fu=fu,ap=ap)
	##view(structure['modules'])
	
	## test call to extract modules
	#input_files=['Ca3Al2Si3O12.cif']
	#structure=extract_module(input_files,bondtable)
	#view(structure['modules'])
	
	#############################################################################
	# startup bits for resuming an old calculation 
	#############################################################################
	
	if restart == True:
		# append to previous output.txt file
		o=open("output.txt",'a')
		#check that the structures folder exisits
		if not os.path.isdir("structures"):
			os.mkdir("structures")
		used_backup=False	
		try:
			os.chdir("restart")
			#print(glob.glob("*"))
			#try reloading files from the restart folder
			
			initial_population=pickle.load(open("initial_population.p",'rb'))#
			using_prebuilt=pickle.load(open("using_prebuilt.p",'rb'))#
			#if prebuilt structures are present, use those
			if os.path.isfile("pre_built_structures.p"):#
				pre_built_structures=pickle.load(open("pre_built_structures.p",'rb'))#
			r=pickle.load(open("r.p",'rb'))#
			ca=pickle.load(open("ca.p",'rb'))#
			T=pickle.load(open("T.p",'rb'))#
			T0=pickle.load(open("T0.p",'rb'))#
			sa=pickle.load(open("sa.p",'rb'))#
			search=pickle.load(open("search.p",'rb'))#
			energies=pickle.load(open("energies.p",'rb'))#
			search_generation_complete=pickle.load(open("search_generation_complete.p",'rb'))#
			
			graph_output=pandas.read_csv("graph_output.csv")#
			graph_output=graph_output.to_dict(orient='list')#
			
			if search == 1 or 2:
				moves=pickle.load(open("moves.p",'rb'))#
				try: #this may not exist yet!
					used_moves=pickle.load(open("used_moves.p",'rb'))
				except:
					used_moves=[]
			#try this, as the search generation may not exisit yet!
			try:
				next_generation=pickle.load(open("next_generation.p",'rb'))
			except:
				pass
			
			#let the user know that the restart files were successfully read
			print("Restart files successfully read")
		
		except: # if the above fails, go to the backup folder and attempt to read those files
			try:
				os.chdir("../")
				os.chdir("backup")
				used_backup=True
				#try reloading files from the backup folder
				
				initial_population=pickle.load(open("initial_population.p",'rb'))
				using_prebuilt=pickle.load(open("using_prebuilt.p",'rb'))
				#if prebuilt structures are present, use those
				if os.path.isfile("pre_built_structures.p"):
					pre_built_structures=pickle.load(open("pre_built_structures.p",'rb'))
				r=pickle.load(open("r.p",'rb'))
				ca=pickle.load(open("ca.p",'rb'))
				T=pickle.load(open("T.p",'rb'))
				T0=pickle.load(open("T0.p",'rb'))
				sa=pickle.load(open("sa.p",'rb'))
				search=pickle.load(open("search.p",'rb'))
				energies=pickle.load(open("energies.p",'rb'))
				search_generation_complete=pickle.load(open("search_generation_complete.p",'rb'))
				
				graph_output=pandas.read_csv("graph_output.csv")
				graph_output=graph_output.to_dict(orient='list')

				
				if search == 1 or 2:
					moves=pickle.load(open("moves.p",'rb'))
					try: #this may not exist yet
						used_moves=pickle.load(open("used_moves.p",'rb'))
					except:
						used_moves=[]
			
				#try this, as the search generation may not exisit yet!
				try:
					next_generation=pickle.load(open("next_generation.p",'rb'))
				except:
					pass

			
			
				# additionally, if we're using the backup, we need to clear the structures folder and resore it from the backup
				os.chdir("../structures")
				cifs=glob.glob("*.cif")
				for i in cifs:
					os.remove(i)
				os.chdir("../backup/structures")
				cifs=glob.glob("*.cif")
				for i in cifs:
					shutil.copy(i,"../../structures/.")
				os.chdir("../")
				
				# also need to restore the solutions files to the main directory:
				os.chdir("backup")
				shutil.copy("cubes.p","../.")
				shutil.copy("hexagonal.p","../.")
				shutil.copy("monoclinic.p","../.")
				shutil.copy("othorhombic.p","../.")
				shutil.copy("tetragonal.p","../.")
				os.chdir("../")
	
				#let the user know that the restart files were successfully read
				print("Could not load restart files, reverted to backup")
			except:
				print("failed to load any restart files, aborting job")
				sys.exit()
				
		# check to see if the calculation had been previously converged
		if r >= rmax:
			print("\n\n***The search routine has already converged. Exiting.*** \n\n")
			sys.exit()
		
		os.chdir("../")
				
	#############################################################################
	# startup bits for a fresh calculation 
	#############################################################################
	
	if restart == False:
		# flag to indicate that there are no outstanding energy calculations to do.
		# count since the current global minimum was identified
		r = 0 
		# count since the last downhill move was made
		ca = 0
		# version of ca count used when determinig whether to switch between search rountines
		sa = 0
		# create new output text file
		o=open("output.txt",'w')
		# create a new dictionary to contain all of the relaxed energies
		energies={}
		#default the search gen to true, set to false when there are structures to optimise:
		search_generation_complete=True
		#if search == 1 or 2, set the moves object to be the initial set.
		if search == 1 or 2:
			moves=n_moves
		#keep track of the moves that we've used:
		used_moves=[]
		
		#set up the object which will make our graph csv file
		graph_output={'move':[],'type':[],'step':[],'energy':[],'temperature':[],'current energy':[],'global minimum energy':[],'structure file name':[],'accepted?':[] }
		
		# create directory for output structures, if it already exisits, clear it out
		if write_all_structures == True:
			if not os.path.isdir("structures"):
				os.mkdir("structures")		
		
		os.chdir("structures")
		if glob.glob("*.cif") != []:
			if platform.system()=='Windows':
				os.system("del *.cif")
			if platform.system()=='Linux':
				os.system("rm *.cif")
		os.chdir("../")
		
		#check to see if there's a restart folder, if so, clear it
		if not os.path.isdir("restart"):
			os.mkdir("restart")
		
		os.chdir("restart")
		to_remove=glob.glob("*")
		if platform.system()=='Windows':
			for i in to_remove:
				os.remove(i)
		if platform.system()=='Linux':
			os.system("rm *")
		os.chdir("../")
		
		#check to see if there is a backup directory, if so clear it
		if not os.path.isdir("backup"):
			os.mkdir("backup")
			
		os.chdir("backup")
		to_remove=glob.glob("*")
		if platform.system()=='Windows':
			for i in to_remove:
				if os.path.isfile(i):
					os.remove(i)
			
			if os.path.isdir("structures"):
				os.chdir("structures")
				to_remove2=glob.glob("*")
				for i in to_remove2:
					os.remove(i)
				os.chdir("../")
			
		if platform.system()=='Linux':
			os.system("rm -r *")
		os.chdir("../")		
		
		# starting here, if it's been selected in the input file, go through and run gn-boss to generate starting structures
		
		#generate_gn_boss_structures=True, # if set to true, when FUSE is firt launched, it will run gn-boss to generate the pool of referennce structures for this calculation.
		#gn_boss_command="C:\ProgramData\anaconda3\_conda.exe run -p C:\ProgramData\anaconda3\envs\ML_FUSE python .\get_cifs_for_FUSE.py", # For my machine, I've setup gn-boss in a seperate python environment, this is the command for that version of python
		#gn_search='tpe', # 'rand' random search, 'tpe' baysian opt, 'pso' particle swarm
		#gn_max_step=5000, # number of generation attempts for gn-boss
		#gn_template_path='C:\Users\cc0u5\Documents\csp_installs\ML-FUSE\template', #path to template files for using gn-boss 
		#gn_zn_range=[1,4], # numbers of formula units to scan with GN-BOSS for generating structures
		#rank_gn_structures=True, #if true, rank the output of the model using SPPs. this is needed if pull_spp_rank = True above		
		#clear_previous_gn_structures=True, #if set to True, before starting the calcluation, remove any previous structures from reference structures & gn-boss generated_results.
		
		if generate_gn_boss_structures == True:
			print("Generating structure pool using Gn-Boss ML model")
			o.write("\nGenerating structure pool using Gn-Boss ML model\n")
			
			
			if clear_previous_gn_structures == True:
				if os.path.isdir("gn-boss"):
					#os.chdir("gn-boss")
					shutil.rmtree("gn-boss") # clear out the previous run.
					#os.chdir("../")
			
					os.chdir(path_to_structures)
					
					f_to_remove=glob.glob("*")
					for w in f_to_remove:
						os.remove(w)
						
					os.chdir("../")
					os.mkdir("gn-boss")
					os.chdir("gn-boss")
			
				else:
					os.mkdir("gn-boss")
					os.chdir("gn-boss")
			
			if clear_previous_gn_structures != True:
				if not os.path.isdir("gn-boss"):
					os.mkdir("gn-boss")
					os.chdir("gn-boss")
			
			#go fetch the template file
			#try:
			shutil.copytree(gn_template_path,'.',dirs_exist_ok=True)
				#os.system("cp -r "+gn_template_path+"\* .")
			#except:
				#os.system("cp -r "+gn_template_path+"* .")
						
			os.chdir("chemical_compositions")
			template=open("template.in",'r').readlines()
			run_files=[]
			#for each value of Z build the input file
			for z in range(gn_zn_range[0],gn_zn_range[-1]+1):
				run_file=template.copy()
				form='compound = '
				total=0
				for w in list(composition.keys()):
					form+=w
					form+=str(composition[w]*(z))
					total+=composition[w]*(z)
					form += " "
				form += "\n"
				if total <= max_atoms:
					run_file[2]=form
				
				run_file[35]="algorithm = "+gn_search+"\n"
				run_file[39]="max_step = "+str(gn_max_step)+"\n"
				
				run_file2=open("gnoa-input_"+str(z)+".in",'w')
				for w in run_file:
					run_file2.write(w)
					
				run_file2.close()
				
				run_files.append("gnoa-input_"+str(z)+".in")
			
			os.remove("template.in")
			
			os.chdir("../")
			#print(os.getcwd())
			#sys.exit()
			#now go and run the calculations
			os.system(gn_boss_command)
			
			#collate the results in the reference structures folder
			os.chdir("results")
			r_files=glob.glob("*")
			for w in r_files:
				if os.path.isdir(w):
					if not w == "best_structures":
						os.chdir(w)
						#print(os.getcwd())
						#sys.exit()
						os.chdir("structures")
						#print(os.getcwd())
						#sys.exit()
						to_copy=glob.glob("*.cif")
						for v in to_copy:
							if v != 'temp.cif':
								shutil.copy(v,"../../../../"+path_to_structures+"/.")
						os.chdir("../../")
						
			os.chdir("../../")
			
			#print(os.getcwd())
			#sys.exit()
			
			#if required, go through and rank the structures
			if rank_gn_structures != None:
				os.chdir(path_to_structures)
				
				r_cifs=glob.glob("*.cif")
				r_results={'file':[],'energy':[],'atoms':[],'converged':[]}
				
				if not os.path.isfile("dummy.lib"):
					fr=open("dummy.lib",'w')
					fr.close()
				
				#get the required spp library
				temp=read(r_cifs[0])
				r_elements=[]
				for w in range(len(temp)):
					if not temp[w].symbol in r_elements:
						r_elements.append(temp[w].symbol)
						
				#print(r_elements)
				#sys.exit()
				assemble_spp(r_elements,spp_path=spp_path)
				
				print("Ranking structures from Gn-Boss ML model using SPPs")
				o.write("\nRanking structures from Gn-Boss ML model using SPPs")
				
				for w in range(len(r_cifs)):
					print(str(w+1)+" of: "+str(len(r_cifs)),end='\r')
					try:
						atoms=read(r_cifs[w])
					except:
						continue
					if ranking == 'gulp':
						try:
							atoms,energy,converged=run_gulp(atoms=atoms,shel=shel,kwds=r_kwds,opts=r_gulp_opts,lib=r_lib,produce_steps=False,gulp_command=gulp_command,gulp_timeout=gulp_timeout)
						except:
							converged=False
							energy=1.e20
					if ranking == 'chgnet':
						try:
							atoms,energy,converged = run_chgnet(atoms,n_opts=n_opts,rel=rel,relaxer_opts=relaxer_opts,opt_class=opt_class,mode=rank_gn_structures,opt_device=opt_device)
						except:
							converged=False
							energy=1.e20

					if ranking == 'mixed':
			
						try:
							atoms,energy,converged=run_calculators(atoms=atoms,vasp_opts=
							vasp_opts,kcut=kcut,produce_steps=None,shel=shel,
							kwds=r_kwds,gulp_opts=r_gulp_opts,lib=r_lib,calcs=r_calcs,dist_cutoff=dist_cutoff,qe_opts=qe_opts,
							gulp_command=gulp_command,gulp_timeout=gulp_timeout,
							n_opts=n_opts,rel=rel,relaxer_opts=relaxer_opts,opt_class=opt_class,
							opt_device=opt_device,mode=rank_gn_structures)
							
						except:
							converged = False
							energy = 1.e20						

					energy=energy/len(atoms)
        
					r_results['file'].append(r_cifs[w])
					r_results['energy'].append(energy)
					r_results['atoms'].append(atoms)
					r_results['converged'].append(converged)
        			
					write(r_cifs[w],atoms)
        			
				dat=pandas.DataFrame.from_dict(r_results).sort_values(['energy'],axis=0,ascending=True)
				dat2=dat.to_dict(orient='list')
				#print(dat2.keys())
				table={'file':[],'energy':[]}
				for j in range(len(dat2[list(dat2.keys())[0]])):
				    table['file'].append(dat2['file'][j])
				    table['energy'].append(dat2['energy'][j])
				    
				dat3=pandas.DataFrame.from_dict(table)
				dat3.to_csv("ranking.csv",index=None)
				    
				os.chdir("../")	
				
			t2=datetime.datetime.now()	
			if generate_structures_only == True:
				print("structure generation complete")
				print("\ntotal time: "+str(t2-t1)+" hours:minutes:seconds")
				o.write("\ntotal time: "+str(t2-t1)+" hours:minutes:seconds\n")
				sys.exit()

			else:
				print("structure generation complete")
				o.write("\ntotal time: "+str(t2-t1)+" hours:minutes:seconds\n")

        #if j == 1:
           #break					
			
			#next bit will be to go through and rank the results!
				
	
	#############################################################################
	# startup bits which we need to do irrespective of restart state
	#############################################################################
	
	# write header to file
	
	o.write("################################################################################")
	o.write("\n#									       #")
	o.write("\n#			   Flexible Unit Structure Engine		       #")
	o.write("\n#				     (FUSE)				       #")
	o.write("\n#				     v2.02				       #")
	o.write("\n#									       #")
	o.write("\n################################################################################")
	o.write("\n\n")

	#############################################################################

	#############################################################################
	# set the system type, for the moment, always forcing the system to neutral #
	#############################################################################
	
	neutral=1
	system_type="neutral"

	#############################################################################

	### write to output file if restarts have been read #########################
	if restart == True:
		if used_backup == False:
			o.write("\nRestart files successfully read")
		if used_backup == True:
			o.write("\nCould not load restart files, reverted to backup")
	#############################################################################

	#### get the composition from the input file to write to output #############
	keys=list(composition.keys())
	string=""
	for i in range(len(keys)):
		string+=keys[i]
		if composition[keys[i]] != 1:
			string+=str(composition[keys[i]])
		
	#############################################################################
	# various write commands from FUSE107 carried over:
	#############################################################################
	
	print ("Input formula = "+string+"\nMaximum number of formula units = "+str(max_fus))
	o.write("\nInput formula = "+string+"\nMaximum number of formula units = "+str(max_fus))
	system_type_string=str( str(system_type) + " input formula" )
	print (system_type_string)
	o.write("\n"+system_type_string)
	print("ap calculated to be: "+str(ap)+" Angstroms")
	o.write("\nap calculated to be: "+str(ap)+" Angstroms")
	
	#### check the search type to be used #######################################
	if swap_searches == True:
		print("Search routine: Mixed")
		o.write("\nSearch routine: Mixed")
		
	if swap_searches == False:
		if search == 1:
			print ("Search routine: Basin Hopping")
			o.write("\nSearch routine: Basin Hopping")
			
		#At the moment, these two haven't been implimented	
		if search == 2:
			print ("Search routine: Basin Hopping with reinforcement learning")
			o.write("\nSearch routine: Basin Hopping with reinforcement learning")
		if search == 3:
			print ("Search routine: Genetic Algorithm")
			o.write("\nSearch routine: Genetic Algorithm")	
	
	#############################################################################
	
	#############################################################################
	# Generate the inintial population of structures
	#############################################################################
	
	# Only generate an initial population if we are starting a fresh calculation
	
	if restart == False:
		generation_complete=False
		print ("\n\n############################ Generating Initial Population ############################\n\n")
		o.write("\n\n############################ Generating Initial Population ############################\n")
		
		initial_population={}
		built=0
		pre_built_structures={}
		# flag to use in the generation loop below to first fill the population with user supplied structures, then randomly generated ones
		using_prebuilt = False
		pre_built_used=1.e20
		# if flagged in input file, go and find the pre-built structures and how many of them there are
		if read_exisiting_structures == True:
			print ("Looking for user supplied structures")
			o.write("\nLooking for user supplied structures")
			os.chdir(path_to_structures)
			
			cifs=glob.glob("*.cif")

			if len(cifs) > 0:
				#if we want to pull structures in a random order do this:
				if pull_random==True:
					# shuffle the order in which the cifs are read in, so that if we have more reference structures than we can use in the initial population, they are drawn in a random order, with the remainder retained for later
					random.shuffle(cifs)
					
				#if we want to pull pre-built structures by an energy ranking from spps.	
				elif pull_spp_rank == True:
					spp_ranks={'file':[],'energy':[]}
					csvs=glob.glob("*.csv")
					for i in csvs:
						temp_dat=pandas.read_csv(i)
						temp_dat_2=temp_dat.to_dict(orient='list')
						for j in range(len(temp_dat_2['file'])):
							spp_ranks['file'].append(temp_dat_2['file'][j])
							spp_ranks['energy'].append(temp_dat_2['energy'][j])
					
					spp_dat=pandas.DataFrame.from_dict(spp_ranks).sort_values(['energy'],axis=0,ascending=True)
					spp_dat2=spp_dat.to_dict(orient='list')
					cifs=spp_dat2['file']					
					
				if len(cifs) == 1:
					print ("... found ",str(len(cifs))," structure")
					o.write("\n... found "+str(len(cifs))+" structure")
				else:
					print ("... found ",str(len(cifs))," structures")
					o.write("\n... found "+str(len(cifs))+" structures")				
				
				pre_built_structures={}
				for i in cifs:
					# flag to determine if we will reject the cif
					reject=False
					# read in the cif
					try:
						temp=read(i)
					except:
						continue
					# first need to check to make sure that the structure does not have partially occupied sites
					occupancies = list(temp.info['occupancy'].values())
					
					for j in occupancies:
						if list(j.values())[0] < 1.0:
							reject=True
							break
					
					# then need to check the composition, does it match the normalised composition from the input file?
					
					symbols=[]
					for j in temp.get_chemical_symbols():
						if not j in symbols:
							symbols.append(j)
					
					symbols.sort()
					input_symbols=list(norm_comp.keys())
					input_symbols.sort()
					
					# first check to see if the list of element symbols matches:
					#print(symbols," ",input_symbols)
					#print(symbols == input_symbols)
					# if it matches, we then need to go through & work out the normalised composition of the structure to see if it matches
					if symbols == input_symbols:
						comp={}
						for j in symbols:
							comp[j]=0
						for j in temp:
							comp[j.symbol]+=1
						
						temp_norm_factor=round( min(list(comp.values())),2)
						
						for j in list(comp.keys()):
							comp[j]=round( comp[j]/temp_norm_factor,2 )
						
						# now go through and compare the formulas, looking for any quantities which do not match
						
						for j in list(norm_comp.keys()):
							if norm_comp[j] != comp[j]:
								reject = True
						
					# if it does not match, then the composition won't match, so can reject the structure at this stage
					if symbols != input_symbols:
						reject = True
					
					# check to see if the size of the read in structure is compatable with the number of atoms in the calculation
					if len(temp) > max_atoms:
						reject = True
					
					# check to see if any distances in the structure are shorter than the dist_cutoff value
					dists=[]
					for j in range(len(temp)):
						boundto=list(range(len(temp)))
						del boundto[j]
						dists=temp.get_distances(j,boundto,mic=True)
						if min(dists) < dist_cutoff:
							reject=True
							break
					
					if reject == True:
						os.remove(i)
					
					#if use_spglib == True:
					#	try:
					#		lattice,positions,numbers=spglib.standardize_cell(temp,symprec=1.e-5)
					#		temp2=Atoms(numbers=numbers,pbc=True)
					#		temp2.cell=lattice
					#		temp2.set_scaled_positions(positions)
					#		temp=temp2.copy()
					#		
					#	except:
					#		pass
					
					#if we haven't found a reason to reject the structure, add it to the pool of pre_build structures, and flag that it hasn't been used yet
					if reject == False:
						pre_built_structures[i]={'structure':temp,'used?':False,'file':i}
				
				if len(list(pre_built_structures.keys())) == 1: 
					print ("... successfully read ",str(len(list(pre_built_structures.keys())))," structure")
					o.write("\n... successfully read "+str(len(list(pre_built_structures.keys())))+" structure")				

				else: 
					print ("... successfully read ",str(len(list(pre_built_structures.keys())))," structures")
					o.write("\n... successfully read "+str(len(list(pre_built_structures.keys())))+" structures")				

				# if we have successfully read in some sample structures, set the flag to let the next part know that we are going to use them!
				if len(list(pre_built_structures.keys())) > 0:
					using_prebuilt=True
					#count that we can use below to keep track of the number of prebuilt structures that we have used
					pre_built_used=0
				
				if len(list(pre_built_structures.keys())) <= 0:
					using_prebuilt=False
					pre_built_used=1.e20
					
			print("")	
			os.chdir("../")
			
		while built < initial_gen:
			#quick check to see if we're still using prebuilt structures:
			if pre_built_used >= len(list(pre_built_structures.keys())):
				using_prebuilt = False
			
			#############################################################################
			# this section for reading in and slicing pre-defined structures from cifs
			#############################################################################
			if using_prebuilt == True:
				if pre_built_used < len(list(pre_built_structures.keys())):
					atoms=pre_built_structures[list(pre_built_structures.keys())[pre_built_used]]['structure']
					#flag that we've now used this structure
					pre_built_structures[list(pre_built_structures.keys())[pre_built_used]]['used?']=True
					#now get the structure in the correct format
					write("temp.cif",atoms)
					structure=extract_module(["temp.cif"],bondtable)
					structure['source']=pre_built_structures[list(pre_built_structures.keys())[pre_built_used]]['file']
				
				pre_built_used +=1
			
			#############################################################################
			# this section for the geneation of random structures 
			#############################################################################
			if using_prebuilt == False:
				structure = get_new_structure(composition=composition,
				max_atoms=max_atoms,imax_atoms=imax_atoms,restart=restart,max_ax=max_ax,
				density_cutoff = density_cutoff,check_bonds= check_bonds,btol= btol,
				check_distances= check_distances,system_type= system_type,
				dist_cutoff = dist_cutoff,vac_ratio = vac_ratio,atoms_per_fu= atoms_per_fu,
				imax_fus= imax_fus,max_fus= max_fus,cubic_solutions= cubic_solutions,
				tetragonal_solutions= tetragonal_solutions,
				hexagonal_solutions= hexagonal_solutions,
				orthorhombic_solutions= orthorhombic_solutions,
				monoclinic_solutions= monoclinic_solutions,bondtable=bondtable,
				ideal_density=ideal_density,fu=fu,ap=ap,use_spglib=use_spglib)
								
			#flag that the structure has not been optimised, and set energy to zero
			structure['optimised?']=False
			structure['energy']=0.0
			structure['converged']=False
			
			#add the structure into the initial generation			
			initial_population[str(built)]=structure
			#print to output a count of the number of structures we have inserted into our initial generation
			print (str("built " + str(built+1) + " of " + str(initial_gen)),end="\r")
			#print(structure)
			#print(structure.keys())
			#sys.exit()
			built +=1
	
	#print an empty line to clear the standard output
	print("                                                                    ", end='\r')
	print("\n")
	o.write("\n")
	#############################################################################
	
	#############################################################################
	# Now need to go through and perform geometry optimisation on the initial 
	# population if any of the structures are not flagged as optimised
	#############################################################################
	
	#check to see if the initial population has been completed, as may have been 
	#in a previous run
	generation_complete=True
	ncomplete=0
	keys_to_complete=[]
	for i in list(initial_population.keys()):
		if initial_population[i]['optimised?']==True:
			ncomplete += 1
			
		else:
			keys_to_complete.append(i)
			
	if ncomplete == len(list(initial_population.keys())):
		generation_complete = True	
	
	if ncomplete < len(list(initial_population.keys())):
		generation_complete = False
	
	#If we've restarted a calculation and continuing the initial population, print this to the output.
	if restart == True:
		if generation_complete == False:
			print("\n\n############################ Resuming Initial Population ############################\n")
			o.write("\n\n############################ Resuming Initial Population ############################\n")
	
	#while we have iterations to do so, go through and continue optimising the initial population
	
		
	while generation_complete == False:
		#first check to see if all are complete
		ncomplete=0
		keys_to_complete=[]
		for i in list(initial_population.keys()):
			if initial_population[i]['optimised?']==True:
				ncomplete += 1
			else:
				keys_to_complete.append(i)
				
		if ncomplete == len(list(initial_population.keys())):
			generation_complete = True
			continue
			
		#before running any structures, create a backup of the restart files
		#print("I am here: ",os.getcwd())
		os.chdir("restart")
		backing=glob.glob("*")
		for i in backing:
			shutil.copy(i,"../backup/"+i)
		os.chdir("../")
		os.chdir("backup")
		if not os.path.isdir("structures"):
			os.mkdir("structures")
		os.chdir("../")
		
		os.chdir("structures")
		backing=glob.glob("*.cif")
		for i in backing:
			shutil.copy(i,"../backup/structures/"+i)
		del backing
		os.chdir("../")
		#make sure we also backup the unit cell shapes files
		backing=glob.glob("*.p")
		for i in backing:
			shutil.copy(i,"backup/.")
		del backing
		
		#print out which structure we're upto
		### change this back to end="\r"
		print (str("structure " + str(ncomplete+1) +" of " + str(len(list(initial_population.keys())))),end="\r")

		# Go through and compute the energies for each of the structures

		#structure we're working on is the first one in the "keys_to_complete" list
		atoms=initial_population[keys_to_complete[0]]['atoms']
		backup_atoms=atoms.copy()
		energy=0.0
		converged=False
		iat=len(atoms)
		t1a=time.time()
		
		converged = False
		energy = 1.e20
		#view(atoms)
		if use_spglib == True:
			try:
				
				lattice,positions,numbers=spglib.standardize_cell(atoms,symprec=1.e-5)
				temp2=Atoms(numbers=numbers,pbc=True)
				temp2.cell=lattice
				temp2.set_scaled_positions(positions)
				atoms=temp2.copy()
				#print("I'm using SPGLIB!")

			except:
				#print("I failed at using SPGLIB!!")
				pass
		#view(atoms)
		#sys.exit()
		if ratt_dist > 0:
			atoms.rattle(ratt_dist)

		
		if ctype == 'gulp':
			try:
				atoms,energy,converged=run_gulp(atoms=atoms,shel=shel,kwds=kwds,opts=gulp_opts,lib=lib,produce_steps=False,gulp_command=gulp_command,gulp_timeout=gulp_timeout)
			except:
				converged = False
				energy = 1.e20
		
		if ctype == 'vasp':
			try:
				atoms,energy,converged=run_vasp(atoms=atoms,vasp_opts=vasp_opts,kcut=kcut,produce_steps=False,dist_cutoff=dist_cutoff)
			except:
				#print('except')
				converged = False
				energy = 1.e20						
		
		if ctype == 'qe':
			try:
				atoms,energy,converged=run_qe(atoms=atoms,qe_opts=qe_opts,kcut=kcut,produce_steps=False)
			except:
				converged=False
				energy=1.e20
		
		if ctype == 'chgnet':
			try:
				atoms,energy,converged = run_chgnet(atoms,n_opts=n_opts,rel=rel,relaxer_opts=relaxer_opts,opt_class=opt_class,opt_device=opt_device,use_spglib=use_spglib)
				
			except:
				converged=False
				energy=1.e20
			
		if ctype == 'mixed':	
			
			try:
				atoms,energy,converged=run_calculators(atoms=atoms,vasp_opts=
				vasp_opts,kcut=kcut,produce_steps=None,shel=shel,
				kwds=kwds,gulp_opts=gulp_opts,lib=lib,calcs=calcs,dist_cutoff=dist_cutoff,qe_opts=qe_opts,
				gulp_command=gulp_command,gulp_timeout=gulp_timeout,
				n_opts=n_opts,rel=rel,relaxer_opts=relaxer_opts,opt_class=opt_class,
				opt_device=opt_device)
				
			except:
				converged = False
				energy = 1.e20
		
		#make sure we have the same number of atoms
		if len(atoms) != iat:
			converged = False
			energy = 1.e20
			
		#check for any unphysical distances:	
		
		temp_atoms=atoms.repeat([2,2,2])
		#distances=min(get_distances(new_atoms=atoms))
		#print (distances)
		temp1=temp_atoms.get_all_distances()
		temp2=[]
		for at in range(len(temp1)):
			for at2 in range(len(temp1[at])):
				if temp1[at][at2] != 0:
					temp2.append(temp1[at][at2])
		distances=min(temp2)
		#print(distances,end=' ')
		if distances <= dist_cutoff:
			converged = False
			energy = 1.e20
		
		
		
		energy = energy/len(atoms)
		energy = float(Decimal(energy).quantize(Decimal(str(e_prec))))
		t2a=time.time()
		initial_population[keys_to_complete[0]]['atoms']=atoms
		
		#update the modules etc. to reflect the optimised atoms object
		write("temp.cif",atoms)
		temp_structure=extract_module(["temp.cif"],bondtable)
		os.remove("temp.cif")
		
		if temp_structure == None: # Usually caused by the calculation going horribly wrong, return the atoms object we started with
			write("temp.cif",backup_atoms)
			temp_structure=extract_module(["temp.cif"],bondtable)
			os.remove("temp.cif")
			converged = False
			energy = 1.e20
		
		
		if revert_atoms == False:	
		   	
			initial_population[keys_to_complete[0]]['modules']=temp_structure['modules']
			initial_population[keys_to_complete[0]]['sub module cell']=temp_structure['sub module cell']
			initial_population[keys_to_complete[0]]['shape in submods']=temp_structure['shape in submods']
			initial_population[keys_to_complete[0]]['nmods']=temp_structure['nmods']
			initial_population[keys_to_complete[0]]['ap']=temp_structure['ap']
			
		initial_population[keys_to_complete[0]]['time']=round(t2a-t1a,1)
		initial_population[keys_to_complete[0]]['energy']=energy
		initial_population[keys_to_complete[0]]['optimised?']=True
		initial_population[keys_to_complete[0]]['converged']=converged
		energies[keys_to_complete[0]]=energy
		
		o.write(str("\n"+str(keys_to_complete[0])+"	 {0:=6.6n}   eV/atom".format(float(energy)).rjust(25)))
		o.flush()
		#write out the structure we have just relaxed
		#try:
		write(str("structures/"+str("I-")+str('{0:0=7n}'.format(int(keys_to_complete[0]))+".cif")),atoms)
		#os.chdir("structures")
		#print(keys_to_complete[0])
		#write("test.cif",atoms)
		#os.chdir("../")
			
		#except:
		#	print("fail")		
		
		#graph_output={'move':[],'type':[],'step':[],'energy':[],'temperature':[],'current energy':[],'global minimum energy':[],'structure file name':[],'accepted?':[] }		
		graph_output['move'].append('initial')
		graph_output['type'].append('I')
		graph_output['step'].append(keys_to_complete[0])
		graph_output['energy'].append(energy)
		graph_output['temperature'].append(T)
		graph_output['current energy'].append(0.0)
		graph_output['global minimum energy'].append(0.0)
		graph_output['structure file name'].append(str( str("I-")+str('{0:0=7n}'.format(int(keys_to_complete[0]))+".cif") ) )
		graph_output['accepted?'].append("-")
		
		pandas.DataFrame.from_dict(graph_output).to_csv("restart/graph_output.csv",index=None)

		
		#write out output / restart files for all of the variables we've got so far:
		os.chdir("restart")
		#things to write out: search, what else? will need to create something to plot out energy graphs as well?
		pickle.dump(initial_population,open("initial_population.p",'wb'))
		pickle.dump(using_prebuilt,open("using_prebuilt.p",'wb'))
		try:
			pickle.dump(pre_built_structures,open("pre_built_structures.p",'wb'))
		except:
			pass
		pickle.dump(r,open("r.p",'wb'))
		pickle.dump(ca,open("ca.p",'wb'))
		pickle.dump(T,open("T.p",'wb'))
		pickle.dump(T0,open("T0.p",'wb'))
		pickle.dump(sa,open("sa.p",'wb'))
		pickle.dump(search,open("search.p",'wb'))
		pickle.dump(energies,open("energies.p",'wb'))
		search_generation_complete=True # make sure that we always go onto start a new generation
		pickle.dump(search_generation_complete,open("search_generation_complete.p",'wb'))
		if search == 1 or 2:
			pickle.dump(moves,open("moves.p",'wb'))

		os.chdir("../")
		
		itr +=1
		
		#break the while loop if we have reached the maximum number of iterations that we want to do
		if itr >= iterations:
			break

	###################################################################################################################
	#	If initial generation complete, print a summary of the lowest energy & which structure it is ###################
	###################################################################################################################
	if generation_complete == True:
		print("Generation completed: ",end=' ')
		o.write("\n\nGeneration completed: ")
		print("Lowest energy structure is: " + str(list(energies.values()).index(min(list(energies.values())))) + " with energy: "+str(round(min(list(energies.values())),6)) +" eV/atom")
		o.write("Lowest energy structure is: " + str(list(energies.values()).index(min(list(energies.values())))) + " with energy: "+str(round(min(list(energies.values())),6)) +" eV/atom\n")
		
		#write out the current best structure as a cif
		write("current_structure.cif",initial_population[str(list(energies.values()).index(min(list(energies.values()))))]['atoms'])
		write("global_minimum.cif",initial_population[str(list(energies.values()).index(min(list(energies.values()))))]['atoms'])
		
		### if search == 1 or 2, set the current structure to be the lowest from the initial generation
		if search == 1 or 2:
			current_structure = initial_population[str(list(energies.values()).index(min(list(energies.values()))))].copy()
		#set the blank next_geneation dictionary
		next_generation={}
		prev_min=min(list(energies.values()))
		calc_converged=False

######################################################################################################################
# Move into using the search routines ################################################################################
######################################################################################################################
##Unlike the previous version of FUSE, we will just continue adding structures to the "initial_population" in order to
#keep everything together

			
	print  ("\n\n###########################  Running search rountine  #################################\n\n")
	o.write("\n\n###########################  Running search rountine  #################################\n\n")
	
	while itr < iterations: # need to use this to start running through the new structures.
				
		if generation_complete == True: # only run the search routine if the initial generation has been completed
			
			#check to see if the previous generation has been built & run
			# if true, then we need to go through and build some new structures, if false, then need to skip to next section ans start optimising structures.
			# once the generation of structures has been created, set search_generation_complete = False to trigger the energy calculators.
			if search_generation_complete == True :
				# Before we do anything else, merge the generation which has just been completed into the main population
				#for x in list(next_generation.keys()):
				#	initial_population[x]=next_generation[x].copy()	
            #
				
				# tell the system that we don't have anything built yet
				built = 0
				# start a new dictionary object which will be formatted as per the main population, to fill with new structures
				next_generation={}
				
				
				
				#first of all, get the structure number that we're upto from the initial_population object
				structure_numbers= list(initial_population.keys())
				for x in range(len(structure_numbers)):
					structure_numbers[x]=int(structure_numbers[x])
					
				#set the structure number to start from for this generation of structures
				structure_number = max(structure_numbers)+1
								
				if search == 1: # Basin hopping routine
					#set about building the new generation
					while built < search_gen_bh: # make move actions until we've built the next geneation
						#when generating a structure, these are the things that need assembling:
						#['modules', 'sub module cell', 'shape in submods', 'nmods', 'ap', 'atoms', 'optimised?', 'energy', 'converged']
						#print(current_structure)
						start_atoms=len(current_structure['atoms'])
						start_fus=start_atoms/sum(list(composition.values()))
						print("start atoms: ",start_atoms)
						print("start fus: ",start_fus)
						
						#next_generation[str(structure_number)],move,pre_built_structures=make_basin_move(
						trial,move,pre_built_structures=make_basin_move(
						
							current_structure,
							moves,
							bondtable,
							grid_spacing,
							exclusion,
							#variables needed for the error checking part of the function
							ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,check_distances,dist_cutoff,
							#bits for changing unit cell shape:
							cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,monoclinic_solutions,
							max_atoms,max_ax,vac_ratio,using_prebuilt,pre_built_structures,max_fus,
							atoms_per_fu,
							imax_atoms,use_spglib,initial_population
							)
						
						print("end atoms", len(trial['atoms']) )
						print("end fus",len(trial['atoms'])/sum(list(composition.values())))
						t_atoms=trial['atoms'].copy()
						trial_fus=len(t_atoms)/sum(list(composition.values()))
						
						counts=[]
						correct=False
						if float(trial_fus).is_integer():
							#now need to check each species
							symbols=t_atoms.get_chemical_symbols()
							for x in list(composition.keys()):
								num1=symbols.count(x)
								counts.append(num1/composition[x])
							if all(x == trial_fus for x in counts):
								correct = True
						
						if correct == True:
							next_generation[str(structure_number)]=trial.copy()
						else:
							continue
						
						
						structure_number+=1
						used_moves.append(move)
						
						print("\nmove :",move,": ",move_des[move])
						o.write(str("\nmove :"+str(move)+": "+str(move_des[move])) )
						built+=1
					
					search_generation_complete = False
					
				if search == 2: # Basin hopping routine with reinforcement learning
					pass
				
				if search == 3: # Genetic Algo
					pass
   		
			if search_generation_complete == False: # while the generation is not completed, go through and find the next structure to optimise & do it!
				# first check to see if all structures are complete
				ncomplete=0
				keys_to_complete=[]
				for x in list(next_generation.keys()):
					if next_generation[x]['optimised?']==True:
						ncomplete+=1
					else:
						keys_to_complete.append(x)
				
				if ncomplete == len(list(next_generation.keys())):
					search_generation_complete = True
					continue
					
				#before running any structures, create a backup of the restart files
				#print("I am here: ",os.getcwd())
				os.chdir("restart")
				backing=glob.glob("*")
				for i in backing:
					shutil.copy(i,"../backup/"+i)
				os.chdir("../")
				os.chdir("backup")
				if not os.path.isdir("structures"):
					os.mkdir("structures")
				os.chdir("../")
				
				os.chdir("structures")
				backing=glob.glob("*.cif")
				for i in backing:
					shutil.copy(i,"../backup/structures/"+i)
				del backing
				os.chdir("../")
				
				#print out which structure we're upto
				### change this back to end="\r"
				#print (str("structure " + str(ncomplete+1) +" of " + str(len(list(next_generation.keys())))),end="\r")
				print (str("structure " + str(ncomplete+1) +" of " + str(len(list(next_generation.keys())))),end="\n")

				atoms=next_generation[keys_to_complete[0]]['atoms']
				backup_atoms=atoms.copy()
				energy=0.0
				converged=False
				iat=len(atoms)
				t1a=time.time()

				if use_spglib == True:
					try:
						lattice,positions,numbers=spglib.standardize_cell(atoms,symprec=1.e-5)
						temp2=Atoms(numbers=numbers,pbc=True)
						temp2.cell=lattice
						temp2.set_scaled_positions(positions)
						atoms=temp2.copy()
						
					except:
						pass

				if ratt_dist > 0:
					atoms.rattle(ratt_dist)


				
				if ctype == 'gulp':
					try:
						print("len atoms start:",len(atoms))
						atoms,energy,converged=run_gulp(atoms=atoms,shel=shel,kwds=kwds,opts=gulp_opts,lib=lib,produce_steps=False,gulp_command=gulp_command,gulp_timeout=gulp_timeout)
						print("len atoms end:",len(atoms))
					except:
						converged = False
						energy = 1.e20
				
				if ctype == 'vasp':
					try:
						atoms,energy,converged=run_vasp(atoms=atoms,vasp_opts=vasp_opts,kcut=kcut,produce_steps=False,dist_cutoff=dist_cutoff)
					except:
						#print('except')
						converged = False
						energy = 1.e20						
				
				if ctype == 'qe':
					try:
						atoms,energy,converged=run_qe(atoms=atoms,qe_opts=qe_opts,kcut=kcut,produce_steps=False)
					except:
						converged=False
						energy=1.e20

				if ctype == 'chgnet':
					try:
						atoms,energy,converged = run_chgnet(atoms,n_opts=n_opts,rel=rel,relaxer_opts=relaxer_opts,opt_class=opt_class,opt_device=opt_device)
						
					except:
						converged=False
						energy=1.e20


				if ctype == 'mixed':	
					
					try:
						atoms,energy,converged=run_calculators(atoms=atoms,vasp_opts=
						vasp_opts,kcut=kcut,produce_steps=None,shel=shel,
						kwds=kwds,gulp_opts=gulp_opts,lib=lib,calcs=calcs,dist_cutoff=dist_cutoff,qe_opts=qe_opts,
						gulp_command=gulp_command,gulp_timeout=gulp_timeout,
						n_opts=n_opts,rel=rel,relaxer_opts=relaxer_opts,opt_class=opt_class,
						opt_device=opt_device)
						
					except:
						converged = False
						energy = 1.e20
				
				#check to make sure we still have the correct number of atoms!
				if len(atoms) != iat:
					converged = False
					energy = 1.e20
				
				#check for any unphysical distances:	
								
				temp_atoms=atoms.repeat([2,2,2])
				#distances=min(get_distances(new_atoms=atoms))
				#print (distances)
				temp1=temp_atoms.get_all_distances()
				temp2=[]
				for at in range(len(temp1)):
					for at2 in range(len(temp1[at])):
						if temp1[at][at2] != 0:
							temp2.append(temp1[at][at2])
				distances=min(temp2)
				#print(distances,end=' ')
				if distances <= dist_cutoff:
					converged = False
					energy = 1.e20

				
				energy = energy/len(atoms)
				energy = float(Decimal(energy).quantize(Decimal(str(e_prec))))	
				t2a=time.time()
				next_generation[keys_to_complete[0]]['atoms']=atoms
				energies[keys_to_complete[0]]=energy
				#update the modules etc. to reflect the optimised atoms object
				
				if revert_atoms==False:
					try:
						write("temp.cif",atoms)
						temp_structure=extract_module(["temp.cif"],bondtable)
						os.remove("temp.cif")
					except:
						temp_structure = None
						
					if temp_structure == None: # Usually caused by the calculation going horribly wrong, return the atoms object we started with
						write("temp.cif",backup_atoms)
						temp_structure=extract_module(["temp.cif"],bondtable)
						os.remove("temp.cif")
						converged = False
						energy = 1.e20
				 
				 
				 
					next_generation[keys_to_complete[0]]['modules']=temp_structure['modules']
					next_generation[keys_to_complete[0]]['sub module cell']=temp_structure['sub module cell']
					next_generation[keys_to_complete[0]]['shape in submods']=temp_structure['shape in submods']
					next_generation[keys_to_complete[0]]['nmods']=temp_structure['nmods']
					next_generation[keys_to_complete[0]]['ap']=temp_structure['ap']
					
				next_generation[keys_to_complete[0]]['time']=round(t2a-t1a,1)
				next_generation[keys_to_complete[0]]['energy']=energy
				next_generation[keys_to_complete[0]]['optimised?']=True
				next_generation[keys_to_complete[0]]['converged']=converged

				#Print to output file where we're up to:
				o.write(str("\n"+str(keys_to_complete[0])+"	 {0:=6.6n}   eV/atom".format(float(energy)).rjust(25)))
				#write out the cif for the structure we've just relaxed
				try:
					write(str("structures/"+str("S-")+str('{0:0=7n}'.format(int(keys_to_complete[0]))+".cif")),atoms)
				except:
					write(str("structures/"+str("S-")+str('{0:0=7n}'.format(int(keys_to_complete[0]))+".cif")),backup_atoms)
				#set this so that we can record it for the graph:
				accepted=False
				#check to see if the generation has been completed & we need to try the acceptance
				
				temp_ncomplete=0
				temp_keys_to_complete=[]
				for x in list(next_generation.keys()):
					if next_generation[x]['optimised?']==True:
						temp_ncomplete += 1
					else:
						temp_keys_to_complete.append(i)
						
				if temp_ncomplete == len(list(next_generation.keys())):
					#first work out the minimum energy:
					lowest=[0,1.e20]
					for x in list(next_generation.keys()):
						if next_generation[x]['energy'] < lowest[1]:
							lowest[1]=next_generation[x]['energy']
							lowest[0]=x
					
					#need to work out the minimum energy prior to this generation:
					tes=[]
					for x in list(initial_population.keys()):
						tes.append(initial_population[x]['energy'])
					prev_min=min(tes)
					
					tes = []
					
					
					gen_min={list(next_generation.keys())[0] : next_generation[list(next_generation.keys())[0]].copy() }
						
					for x in list(next_generation.keys()):
						if next_generation[x]['energy'] < gen_min[list(gen_min.keys())[0]]['energy']:
							gen_min ={x: next_generation[x].copy()}
					
					#print ("generation minimum: ",gen_min)
					
					
					del tes
					dE_glob=gen_min[list(gen_min.keys())[0]]['energy']-prev_min
					dE_glob=float(Decimal(str(dE_glob)).quantize(Decimal(str(e_prec))))
					
					dE_curr=gen_min[list(gen_min.keys())[0]]['energy']-current_structure['energy']
					dE_curr=float(Decimal(str(dE_curr)).quantize(Decimal(str(dE_curr))))
				
					# MC accept / reject step:
					
					accept_move=False
					accepted=False
					
					if dE_curr < 0.: # if we have a downhill move, we have to accept it
						accept_move=True
						accepted=True
						current_structure=gen_min[list(gen_min.keys())[0]].copy()
						write("current_structure.cif",current_structure['atoms'])
		

						#reset the counts since a downhill move
						ca=0
						sa=0
						
						#reset the temperature parameter
						T=T0
						
						#reset moves
						moves=n_moves
						
						#if it's the new global minimum, reset the r count
						if dE_glob < 0.:
							r=0
							write("global_minimum.cif",current_structure['atoms'])
					
					#only test this part if the energy is higher than current
					if dE_curr >= 0.:
						if dE_curr > 0.:
							rand=random.random()
							Test=math.exp(-dE_curr/T)
							#If this, accept the new structure
							if Test >= rand:
								accept_move=True
								accepted=True
								current_structure=gen_min[list(gen_min.keys())[0]].copy()
								write("current_structure.cif",current_structure['atoms'])
						
						#update r, ca and sa
						r+=search_gen_bh
						ca+=search_gen_bh
						sa+=search_gen_bh
					
					#check to see if the system needs to be melted
					if accept_move == False:
						if ca >= melt_threshold:
							T += random.random()/100
							moves=r_moves
					
					#check to see if the BH rountine has comverged:
					calc_converged=False
					if r >= rmax:
						calc_converged=True
					
					#print(current_structure)
					
					
					
					
					print("Generation completed")
					o.write("\nGeneration completed")
					print(str("E = "+str("{0: .5e}").format(float(lowest[1])) + " dE vs. global = " + str("{0: .4e}").format(float(dE_glob)).rjust(7)+" r: "+str(r).rjust(4))+"	 T = "+str("{0:.5f}").format(T),end='\n')
					o.write(str("\nE = "+str("{0: .5e}").format(float(lowest[1])) + " dE vs. global = " + str("{0: .4e}").format(float(dE_glob)).rjust(7)+" r: "+str(r).rjust(4))+"	 T = "+str("{0:.5f}").format(T)+"\n")
					o.flush()
					
					# insert the structures from this generation into the main population
					for x in list(next_generation.keys()):
						initial_population[x]=next_generation[x].copy()	
				
				
				#graph_output={'move':[],'type':[],'step':[],'energy':[],'temperature':[],'current energy':[],'global minimum energy':[],'structure file name':[],'accepted?':[] }		
				graph_output['move'].append(move)
				graph_output['type'].append('S')
				graph_output['step'].append(keys_to_complete[0])
				graph_output['energy'].append(energy)
				graph_output['temperature'].append(T)
				graph_output['current energy'].append(current_structure['energy'])
				graph_output['global minimum energy'].append(prev_min)
				graph_output['structure file name'].append(str( str("S-")+str('{0:0=7n}'.format(int(keys_to_complete[0]))+".cif") ) )
				graph_output['accepted?'].append(str(accepted))

				
				
				#write out output / restart files for all of the variables we've got so far:
				os.chdir("restart")
				#things to write out: search, what else? will need to create something to plot out energy graphs as well?
				pickle.dump(initial_population,open("initial_population.p",'wb'))
				pickle.dump(using_prebuilt,open("using_prebuilt.p",'wb'))
				try:
					pickle.dump(pre_built_structures,open("pre_built_structures.p",'wb'))
				except:
					pass
				pickle.dump(r,open("r.p",'wb'))
				pickle.dump(ca,open("ca.p",'wb'))
				pickle.dump(T,open("T.p",'wb'))
				pickle.dump(T0,open("T0.p",'wb'))
				pickle.dump(sa,open("sa.p",'wb'))
				pickle.dump(search,open("search.p",'wb'))
				pickle.dump(energies,open("energies.p",'wb'))
				pickle.dump(search_generation_complete,open("search_generation_complete.p",'wb'))
				if search == 1 or 2:
					pickle.dump(moves,open("moves.p",'wb'))
					pickle.dump(used_moves,open("used_moves.p",'wb'))
				pickle.dump(next_generation,open("next_generation.p",'wb'))
      		
				pickle.dump(graph_output,open("graph_output.p",'wb'))
				pandas.DataFrame.from_dict(graph_output).to_csv("graph_output.csv",index=None)
      		
				os.chdir("../")
								
				if calc_converged==True:
					print("\n\n*** BH rountine converged, stopping calculation ***\n")
					o.write("\n\n\n*** BH rountine converged, stopping calculation ***\n")
					break
				
				#check to see if a stop file has been written:
				if os.path.isfile("stop.txt"):
					print("\n\n ***** STOP file found, exiting job ***** \n")
					o.write("\n\n\n ***** STOP file found, exiting job ***** \n\n")
					os.remove("stop.txt")
					break
				
				itr+=1





#wrap up things: printout structure number and energy of current minimum & time taken
	print("\n\n################################# Calculation finished ##############################\n")
	o.write("\n\n################################# Calculation finished ##############################\n")
	
	print("Lowest energy structure is: " + str(list(energies.values()).index(min(list(energies.values())))) + " with energy: "+str(round(min(list(energies.values())),6)) +" eV/atom")
	o.write("\nLowest energy structure is: " + str(list(energies.values()).index(min(list(energies.values())))) + " with energy: "+str(round(min(list(energies.values())),6)) +" eV/atom\n")

	#copy the graph output file to the main directory
	try:
		shutil.copy("restart/graph_output.csv",".")
		if output_graph_at_end == True:
			plot_graph(search=1)
	except:
		pass
	
	times=[]
	for x in list(initial_population.keys()):
		if 'time' in list(initial_population[x].keys()):
			times.append(initial_population[x]['time'])
	
	if len(times) > 0:
		import statistics
		mtime=round(statistics.mean(times),1)
		try:
			stdevtime=round(statistics.stdev(times),1)
		except:
			stdevtime=0
	
		print("\nMean per structure: "+str(mtime)+"("+str(stdevtime)+") seconds")
		o.write("\nMean per structure: "+str(mtime)+"("+str(stdevtime)+") seconds\n")
	
	t2=datetime.datetime.now()
	print("\ntotal time: "+str(t2-t1)+" hours:minutes:seconds")
	o.write("\ntotal time: "+str(t2-t1)+" hours:minutes:seconds\n")
	
	#close the output file and exit
	o.close()
	sys.exit()
	
# next things to do:

#Rebuild the search routine(s):
#  come up with the list of actions which need building:
#  list of actions to build for the basin hopping (perhaps think about sorting the order out too:
#	1. swap two atoms - DONE
#	2. swap an atom in to a vacent space - DONE

#	3. swapping the position of more than two atoms - no vacancies - DONE
#	4. swapping the position of more than two atoms - including vacancies - DONE
#	5. swap all of the atoms - no vacancies -DONE
#	6. swap all of the atoms - including vacancies -DONE

#  Think that for both 7 & 8, need to think about lattice translations in the modules in order to avoid atom crashes, at the moment, we seem to be getting alot of atom crashes. other option, when creating the modules, try to flatten them to get them back to where we were with fuse originally?
#	7. swap the positions of two sub-modules - DONE
#	8. find two full slices/modules in the structure and switch their positions - DONE

#	9. generate a new set of instructions for the current modules set; equiv here may be trying to work out different possible shapes, perhaps using the solutions function from fuse107 - DONE

#	10. move inspired by genetic algorithms, e.g. mutating / mating the structure, here may not have to be a full cut between structures, may be able to locat sub-modules or modules that we can swap in - DONE

#	11. attmpt to grow the structure along an axis. here there's also the oppertunity to apply symmetry operations to the grown part of the cell, e.g. a mirror plane or inversion centre - DONE
#  12. triple the length of a structure, akin to move 11, again with the oppertunity to apply symmetry operations to the new part, e.g. translating it by 1/3. each time - DONE

#	13. random new structure with upto the same number of fus as we have currently - DONE
#	14. random new structure, allowing the use of any remaining pre-built structures - DONE


#	basin hopping with RLCSP
#	GA
#	hybrid, using BH and GA
#
#other bits:
#Re instate the graph outputs & csv files for data analysis - DONE with using basic BH routine		
#Tweak the module extraction routine, see what happens if each the co-ordinatese for each sub-module are flattened slightly, e.g. make it so that the maximum value along z = 0.5 or 0.75?