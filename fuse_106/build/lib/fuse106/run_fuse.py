import datetime
import numpy
import sys
from ase import *
from fuse106.create_random_string import *
from fuse106.error_check_structure import *
from fuse106.assemble_structure import *
from fuse106.create_random_instructions import *
from fuse106.possible_solutions import *
from fuse106.run_gulp import *
from fuse106.run_vasp import *
from fuse106.run_qe import *
from ase.visualize import *
import warnings
from decimal import *
import os
import platform
import glob
from fuse106.get_distances import *
import math
from fuse106.make_basin_move import *
import random
import shutil
from fuse106.plot_results import plot_graph
from ase.build import sort
from fuse106.run_multiple_calculators import run_calculators
from fuse106.generate_random_structure import generate_random_structure
import pandas

#-------------------------------------------------------------------------------
#Known issues / to do list:
# - Currently the neutral input formula type does not work with the GA, by the looks
# of it, I never built the function for it! Need to mimic how it is done for ionic
# for neutral (removing the parts which split the structures into anions and cations
#
#
#-------------------------------------------------------------------------------

warnings.filterwarnings("ignore") # currently use this as python raises RuntimeError 
# warnings when a physically unreasonable unit cell is generated, this is caught
# and structures rejected by the "converged" variable, so the warnings aren't really
# needed. Comment this out for debugging.

t1=datetime.datetime.now()

print("################################################################################")
print("#									      #")
print("#		       Flexible Unit Structure Engine			      #")
print("#				 (FUSE)					      #")
print("#				 v1.06            			      #")
print("#				                          		      #")
print("################################################################################")
print("\n\n")



def run_fuse(composition='',search=1,initial_gen='',max_atoms='',vac_ratio=4,ap_scale=1.0,
	 max_ax=40,restart=False,ctype='',kwds='',gulp_opts='',lib='',shel='',vasp_opts='',
	 kcut=30,serial=True,write_all_structures=True,ratt_dist=0.05,density_cutoff=0.4,
	 check_bonds=True,btol=0.25,check_distances=True,dist_cutoff=1.2,iterations=0,
	 search_gen=1,n_moves={1:1,2:1,3:1,4:1,5:1,6:1,7:1,9:1,10:1,11:1},r_moves={3:1,4:1,5:1,6:2,7:3},rmax='',T=0.02,
	 imax_atoms='',produce_steps=False,melt_threshold=500,
	 write_graph=False,
	 new_structure_attempts=5000,calcs='',qe_opts='',gulp_command='gulp < gulp.gin > gulp.got',record_data_log=True,
	 invert_charge=False,write_graph_end=True,gulp_timeout=''):

	
	### run possible solutions to get library of unit cell dimensions in sub-modules
	cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,monoclinic_solutions=possible_solutions(max_ax,restart)
	#############################################################################
	
	### try to read in restart_files ############################################
	T_0=T
	if restart == True:
		try:
			initial_population=pickle.load(open("initial_structures.p",'rb'))
			ini_energies=pickle.load(open("ini_energies.p",'rb'))
			r=pickle.load(open("r.p",'rb'))
			ca=pickle.load(open("ca.p",'rb'))
			if ca >= melt_threshold:
				T=pickle.load(open("T.p",'rb'))
			#if write_graph==True:
			temp=pandas.read_csv("graph_output.csv")
			graph=temp.to_dict(orient='list')
			try:
				del graph['file_name']
			except:
				pass
		
			search=pickle.load(open("search.p",'rb'))
		except:
			restart = False
		 
		#When constructed, will also need to be able to read in results from a search
		#routine
		try:
			search_structures=pickle.load(open("search_structures.p",'rb'))
		except:
			pass #put pass in, as search may not yet have begun, so search structures
		#may not yet exist!
		
		try:
			moves=pickle.load(open("moves.p",'rb'))
		except:
			pass #put pass in, as search may not yet have begun, so search structures
		#may not yet exist!
		
		try:
			generation_moves=pickle.load(open("generation_moves.p",'rb'))
		except:
			generation_moves=[]
		
		try:
			energies=pickle.load(open("energies.p",'rb'))
			ini_energies=pickle.load(open("ini_energies.p",'rb'))
		except:
			restart = False
		
		if len(ini_energies) >= initial_gen:
			initial_complete=1
		else:
			initial_complete=0
			ncalc=len(ini_energies)
			
		if initial_complete==1:
			try:
				current_structure=pickle.load(open("current_structure.p",'rb'))
			except:
				restart=False

		try:
			generation_complete=pickle.load(open("generation_complete.p",'rb'))
		except:
			restart=False
		
		if generation_complete == False:
			try:
				generation=pickle.load(open("generation.p",'rb'))
				generation_energies=pickle.load(open("generation_energies.p",'rb'))
			except:
				restart=False
		
		if initial_complete==1:
			try:
				search_structures=pickle.load(open("search_structures.p",'rb'))
				count=len(list(search_structures.keys()))
			except:
				restart=False
			
			
		if restart == True:
			print ("Restart files read... ")
		else:
			print("Error in reading restart files, starting fresh calculation")
	#############################################################################
	
	### startup bits ############################################################
	if restart == False:		   
		prev_start_points=[]
		ini_energies=[]
		energies=[]
		generation_complete=True
		r=0
		ca=0
		sa=0
		#if write_graph==True:
		graph={}
		T_0=T
		moves=n_moves
	itr = 0 # number of structures computed in current run of code
	
	if produce_steps==True:
		if not os.path.isdir("steps"):
			os.mkdir("steps")
		if restart==False:
			if platform.system() == 'Windows':
				os.chdir("steps")
				files=glob.glob("*.cif")
				for z in range(len(files)):
					os.remove(files[z])
				os.chdir("../")
			if platform.system() == 'Linux':
				try:
					os.system("rm steps/*")
				except:
					pass

	#############################################################################

	#### create output files/folders to write to ################################
	if restart == False:
		o=open("output.txt",'w')
	if restart == True:
		o=open("output.txt",'a')
	
	if write_all_structures == True:	
		if not os.path.isdir("structures"):
			os.mkdir("structures")
		os.chdir("structures")
		
		if restart == False:
			if glob.glob("*.cif") != []:
				if platform.system()=='Windows':
					os.system("del *.cif")
				if platform.system()=='Linux':
					os.system("rm *.cif")	
		os.chdir("../")
		
	o.write("################################################################################")
	o.write("\n#										 #")
	o.write("\n#			   Flexible Unit Structure Engine			 #")
	o.write("\n#				     (FUSE)					 #")
	o.write("\n#				     v1.06					 #")
	o.write("\n#										 #")
	o.write("\n################################################################################")
	o.write("\n\n")
	
	### write to output file if restarts have been read #########################
	if restart == True:
		o.write("\nrestart files read... ")
	#############################################################################
	
	#############################################################################
	
	#### Print out the composition from the input file ##########################
	keys=list(composition.keys())
	string=""
	for i in range(len(keys)):
		string+=keys[i]
		if len(composition[keys[0]])==1:
			if composition[keys[i]][0] != 1:
				string+=str(composition[keys[i]][0])
		
		if len(composition[keys[0]])==2:
			if composition[keys[i]][0] != 1:
				string+=str(composition[keys[i]][0])
	#############################################################################

	####compute maximum number of formula units allowed##########################
	### work out how many atoms per formular unit ###############################
	atoms_per_fu = 0
	for i in range(len(keys)):
		atoms_per_fu += composition[keys[i]][0]
	
	
	max_fus=int(max_atoms/atoms_per_fu)
	imax_fus=int(imax_atoms/atoms_per_fu)
	if max_fus == 0:
		print ("ERROR: maximum number of Fus equals 0!, please increase maximum number of atoms in input!")
		sys.exit()
	
	if imax_fus == 0:
		imax_fus=1		
	
	print ("Input formula = "+string+"\nMaximum number of formula units = "+str(max_fus))
	o.write("\nInput formula = "+string+"\nMaximum number of formula units = "+str(max_fus))
	#############################################################################
	
	### create a string containing the atomic numbers for 1 FU ##################
	fu=[]
	for i in range(len(keys)):
		for j in range(composition[keys[i]][0]):
			fu.append(Atoms(keys[i]).get_atomic_numbers()[0])
	#############################################################################
	
	#### check if the composition is presented as ionic #########################
	if len(composition[keys[0]]) != 1:
		print ("Ionic input formula")
		o.write("\nIonic input formula")
		charge=0
		system_type="ionic"		
		for i in range(len(keys)):
			for j in range(composition[keys[i]][0]):
				charge += composition[keys[i]][1]
				
		if charge == 0:
			neutral=1
			
		if charge != 0:
			neutral=0
			
	if len(composition[keys[0]]) == 1:
		print ("Neutral input formula")
		o.write("\nNeutral input formula")
		neutral=1
		system_type="neutral"
		
	if neutral != 1:
		print ("Error: Non-charge neutral formular input")
		o.write("\nError: Non-charge neutral formular input")
		sys.exit()
	if neutral ==1:
		pass
	#############################################################################
	# try and load bond table depending on which type of system we're in
	if restart == True:
		try:
			bondtable=numpy.load("bondtable.npz",allow_pickle=True)
			bondtable=bondtable['bond_table'].item()
		except:
			if system_type == 'ionic':
				if invert_charge == False:
					import fuse106.bond_table_ionic
					bondtable=numpy.load("bondtable.npz",allow_pickle=True)
					bondtable=bondtable['bond_table'].item()
			
				if invert_charge == True:
					import fuse106.bond_table_ionic_invert
					bondtable=numpy.load("bondtable.npz",allow_pickle=True)
					bondtable=bondtable['bond_table'].item()
				
			
			if system_type == 'neutral':
				import fuse106.bond_table_atomic
				bondtable=numpy.load("bondtable.npz",allow_pickle=True)
				bondtable=bondtable['bond_table'].item()
	if restart == False:
		if system_type == 'ionic':
			if invert_charge == False:
				import fuse106.bond_table_ionic
				bondtable=numpy.load("bondtable.npz",allow_pickle=True)
				bondtable=bondtable['bond_table'].item()
			
			if invert_charge == True:
				import fuse106.bond_table_ionic_invert
				bondtable=numpy.load("bondtable.npz",allow_pickle=True)
				bondtable=bondtable['bond_table'].item()
					
		
		if system_type == 'neutral':
			import fuse106.bond_table_atomic
			bondtable=numpy.load("bondtable.npz",allow_pickle=True)
			bondtable=bondtable['bond_table'].item()
		
	
	#### compute the ap to be used for the sub-modules ##########################
	# convert the fu back to symbols
	symbol_fu=[]
	for i in range(len(fu)):
		temp=Atoms(numbers=[fu[i]])
		symbol_fu.append(temp.get_chemical_symbols()[0])
	# read in shannon radi table	
	temp=numpy.load("bondtable.npz",allow_pickle=True)
	bond_table=temp['bond_table'].item()
	ap=0
	if system_type=="neutral":
		for i in range(len(symbol_fu)):
			average = 0
			temp=list(bond_table[symbol_fu[i]].keys())
			for j in range(len(temp)):
				average+=bond_table[symbol_fu[i]][temp[j]][-1]
			average=average/len(temp)
			ap+=average

	ap=ap/len(symbol_fu)
	ap=4*ap
	
	if system_type=="ionic":
		cat_ap=0
		an_ap=0
		cats=0
		anis=0
		for i in range(len(symbol_fu)):
			sym=symbol_fu[i]
			try:
				if composition[sym][-1] > 0:
					cat_ap+=bond_table[sym][composition[sym][-1]][-1]
					cats+=1
				if composition[sym][-1] < 0:
					an_ap+=bond_table[sym][composition[sym][-1]][-1]
					anis+=1
					
			except KeyError:
				print("***** WARNING: dectected ionic species not included in bond table, treating species as an average of it's bondtable entry *****")
				average=0
				temp=list(bond_table[symbol_fu[i]].keys())
				for j in range(len(temp)):
					average+=bond_table[symbol_fu[i]][temp[j]][-1]
				average=average/len(temp)
				ap+=average
					
		cat_ap = (cat_ap / cats)*2
		an_ap  = (an_ap  / anis)*2
		ap = cat_ap + an_ap

	if ap_scale != '':
		ap=ap*ap_scale
	ap = float(Decimal(ap).quantize(Decimal('1e-4')))
	print("ap calculated to be: "+str(ap)+" Angstroms")
	o.write("\nap calculated to be: "+str(ap)+" Angstroms")
	#############################################################################

	#try to work out maximum density ############################################
	# currently this won't work for neutral compounds ,especially if one or more elements
	# can present as either an anion or cation
	mass=0
	volume=0	
	if system_type == 'neutral':
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
		#print(ideal_density)
	if system_type == 'ionic':
		for i in range(len(fu)):
			try:
				temp=Atoms(numbers=[fu[i]])
				mass+=temp.get_masses()[0]
				temp=temp.get_chemical_symbols()[0]
				sym=temp
				temp=bondtable[temp]
				chg_state=composition[sym][1]
				rs=temp[chg_state][-1]
				volume+=((4/3)*math.pi*(rs**3))
			except KeyError:
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

	#### check the search type to be used #######################################
		
	if search == 1:
		print ("Search routine: Basin Hopping")
		o.write("\nSearch routine: Basin Hopping")
		

	#############################################################################

	#### work on generation of initial population ###############################
	if restart == False:
		r=0
		print ("\n\n############################ Generating Initial Population ############################\n\n")
		o.write("\n\n############################ Generating Initial Population ############################\n")
		initial_population={}
		built=0
		sizes={}
		x=0
		#max_fus=16
		while built < initial_gen:
			#try:
			#max_fus=int(max_atoms/atoms_per_fu)
			#imax_fus=int(imax_atoms/atoms_per_fu)
			#print("\natoms per fu",atoms_per_fu)
			target_fu=choice(list(range(1,imax_fus+1)))
			target_atoms=atoms_per_fu*target_fu
			#print("\nimax_fus: ",imax_fus)
			#print("target atoms:",target_atoms,end="\n")
			len_atoms=0
			#while len_atoms != target_atoms:
				#string,instructions=create_random_string(cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,monoclinic_solutions,atoms_per_fu,fu,vac_ratio=vac_ratio,max_fus=imax_fus,system_type=system_type,composition=composition,ap=ap)
				#instructions=create_random_instructions(string,ap,cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,monoclinic_solutions,instructions)
				#atoms,instructions=assemble_structure(string,instructions)
				#accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=imax_atoms)
				#if accept == 0:
				#	continue
				#if len(atoms) > imax_atoms:
				#	continue	
				#if len(atoms) % len(fu) != 0:
				#	continue
			atoms,string,instructions,accept=generate_random_structure(target_atoms,cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,
				monoclinic_solutions,atoms_per_fu,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,
				vac_ratio=vac_ratio,max_fus=imax_fus,
				target_number_atoms=imax_atoms)
			if len(atoms) == 0:
				print("\n\nZero atoms!!!\n\n")
			len_atoms=len(atoms)
			#print("generated_atoms: ",len_atoms,"\n")
				
			temp_atoms=atoms.copy()
			sort(temp_atoms)
			#print(prev_start_points)
			atoms.rattle(ratt_dist)
			initial_population[built+1]=[atoms,string,instructions]
			print (str("built " + str(built+1) + " of " + str(initial_gen)),end="\r")
			### for debugging, just have a routine here that creates a dictionary to plot
			# the distribution of number of atoms per structure #########################
			try:
				sizes[len(atoms)]+=1
			except:
				sizes[len(atoms)]=1		
			built += 1
			#except:
			#	pass
		#view(atoms)
		#print(sizes)
		#sys.exit()
	#############################################################################
	
	####	now go through and actually do something with the structures!! ########
	if restart == False:
		ncalc=0 # keep track of which structure from the initial population has been computed
		initial_complete = 0 # keep track of if the initialisation loop has been completed	
		search_structures={}
		ini_energies=[]
	if restart == True:
		if initial_complete == 0:
			count=1
			print("\n\n############################ Resuming Initial Population ############################\n")
			o.write("\n\n############################ Resuming Initial Population ############################\n")
	if initial_complete == 0:
		if serial == True: # run the structures one at a time
			#itr=0 # number of structures computed in this loop
			for i in range(0+(ncalc), initial_gen):
				if itr == iterations:
					break			
				nstruct=len(list(initial_population.keys()))
				o.write("\n")
				atoms=initial_population[i+1][0]
				#print(initial_population[i+1])
				iat=len(atoms)
				print (str("structure " + str(i+1) +" of " + str(nstruct)),end="\r")
				
				if ctype == 'gulp':
					try:
						atoms,energy,converged=run_gulp(atoms=atoms,shel=shel,kwds=kwds,opts=gulp_opts,lib=lib,produce_steps=produce_steps,gulp_command=gulp_command,gulp_timeout=gulp_timeout)
					except:
						converged = False
						energy = 1.e20
				
				if ctype == 'vasp':
					try:
						atoms,energy,converged=run_vasp(atoms=atoms,vasp_opts=vasp_opts,kcut=kcut,produce_steps=produce_steps,dist_cutoff=dist_cutoff)
					except:
						#print('except')
						converged = False
						energy = 1.e20						
				
				if ctype == 'qe':
					try:
						atoms,energy,converged=run_qe(atoms=atoms,qe_opts=qe_opts,kcut=kcut,produce_steps=produce_steps)
					except:
						converged=False
						energy=1.e20
						
				if ctype == 'mixed':	
					
					try:
						atoms,energy,converged=run_calculators(atoms=atoms,vasp_opts=
						vasp_opts,kcut=kcut,produce_steps=produce_steps,shel=shel,
						kwds=kwds,gulp_opts=gulp_opts,lib=lib,calcs=calcs,dist_cutoff=dist_cutoff,qe_opts=qe_opts,
						gulp_command=gulp_command,gulp_timeout=gulp_timeout)
						
					except:
						converged = False
						energy = 1.e20
						
				if len(atoms) != iat:
					converged = False
					energy = 1.e20
				if converged != False:
					if produce_steps==True:
						targets=glob.glob("atoms*.cif")
						targets.sort()
						for z in range(len(targets)):
							label=str("steps/I"+str(i+1)+"_"+str(z+1)+str(".cif"))
							shutil.copy(targets[z],label)
				energy = energy/len(atoms)
				energy = float(Decimal(energy).quantize(Decimal('1e-6')))
				initial_population[i+1].append(energy)
				initial_population[i+1].append("ini")								
				initial_population[i+1][0]=atoms
				#print(str("I{0:=6n}".format(i+1)+"  {0:0=6.6n}	  eV/atom".format(float(energy)).rjust(25)))
				o.write(str("I{0:=7n}".format(i+1)+"  {0:=6.6n}	  eV/atom".format(float(energy)).rjust(25)))
				if write_all_structures == True:
					try:
						write(str("structures/"+str("I-")+str('{0:0=7n}'.format(i+1)+".cif")),atoms)
					except:
						pass
				
				#if write_graph==True:
				#if converged != False:
				if i==0:
					graph['move']=['ini']
					graph['type']=["I"]
					graph['step']=[i+1]
					graph['energies']=[energy]
					if search == 1:
						graph['temp']=[T]
					graph['current_energy'] = [0]
				else:
					try:
						graph['move'].append('ini')
						graph['type'].append("I")
						graph['step'].append(i+1)
						graph['energies'].append(energy)
						graph['current_energy'].append(0)
						if search == 1:
							graph['temp'].append(T)
					except:
						graph['move']=['ini']
						graph['type']=["I"]
						graph['step']=[i+1]
						graph['energies']=[energy]
						if search == 1:
							graph['temp']=[T]
						else:
							pass
								#graph['temp']=[T]
						graph['current_energy'] = [0]
													
				graph_to_write=pandas.DataFrame(graph)
				graph_to_write.to_csv(path_or_buf="graph_output.csv",index=False)
				if i % 10 == 0:
					if write_graph==True:
						plot_graph(search=search)

				ini_energies.append(energy)
				energies=[min(ini_energies)]
				if os.path.isfile("stop.txt"):
					print ("** Stopping caclulation **")
					os.remove("stop.txt")
					break
				itr+=1
				#print(initial_population[i+1])
				o.flush()
		if itr + ncalc == initial_gen:
			initial_complete = 1
			#if write_graph==True:
			#	plot_graph()
		if initial_complete == 1:
			energies=[]
			print("Initial population completed")
			o.write(str("\nInitial population completed"))
			print(str("lowest energy initial structure: " + str(ini_energies.index(min(ini_energies))+1) + " , " + str(min(ini_energies)))+" eV/atom")
			o.write(str("\n\nLowest energy initial structure: " + str(ini_energies.index(min(ini_energies))+1) + " , " + str(min(ini_energies)))+" eV/atom")
			current_structure=initial_population[ini_energies.index(min(ini_energies))+1]
			search_structures[0]=initial_population[ini_energies.index(min(ini_energies))+1]
			energies=[search_structures[0][3]]
#			write(str("structures/"+str("BH-")+str('{0:0=7n}'.format(0)+".cif")),atoms)
			write("lowest_energy_structure.cif",search_structures[0][0])
			write("current_structure.cif",search_structures[0][0])
	#############################################################################
	o.close()
	o=open("output.txt",'a')
	### now we can move onto using some of the searching routines ###############
	start=1
	if initial_complete == 1:
		if search != 3:
			while r < rmax:
				if start == 1:
					if restart == False:
						#plot_graph(search)
						count=1
					if search == 1:
						print ("\n\n################################# Basin Hopping search ##################################\n")
						print ("Current lowest energy structure: " + str(energies.index(min(energies))) + " , " + str(min(energies))+" eV/atom\n")
						o.write("\n\n################################# Basin Hopping search ##################################\n")
						o.write("\n\nCurrent lowest energy structure: " + str(energies.index(min(energies))) + " , " + str(min(energies))+" eV/atom\n")
						start=0
											
			
				# first need to generate structures for the following generation ############
				if generation_complete==True:
						x=0
						built = 0
						generation=[]
						generation_energies=[]
						generation_moves=[]
						while built < search_gen:
							#print("built",built)
							if search == 1:
								new_structure=make_basin_move(current_structure,atoms_per_fu,fu,vac_ratio,max_fus,system_type,composition,ap,
									cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,monoclinic_solutions,
									moves=moves,max_atoms=max_atoms,initial_population=initial_population,search_structures=search_structures,ideal_density=ideal_density,
									density_cutoff=density_cutoff,check_bonds=check_bonds,btol=btol,bondtable=bondtable,
									check_distances=check_distances,dist_cutoff=dist_cutoff,imax_atoms=imax_atoms)
								
								move = str(new_structure[-1])
								del new_structure[-1]
								
							#print("\n",len(new_structure),"\n")
							atoms=new_structure[0]
							string=new_structure[1]
							instructions=new_structure[2]
							#move="place_holder"
							
							#sys.exit()
							
							if atoms == None:
								continue
							
							#print("\n\n", str(move), "\n\n")
																												
							#initial_population[built+1]=[atoms,string,instructions]
							print (str("built " + str(built+1) + " of " + str(search_gen)),end="\r")
							temp_atoms=atoms.copy()
							sort(temp_atoms)
							atoms.rattle(ratt_dist)
							generation.append(new_structure)
							if search == 1:
								generation_moves.append(str(move))
							else:
								generation_moves.append("G")
							built += 1
						generation_complete=False
			#############################################################################	
				## now need to go through & relax structures ################################
				#print(generation_moves)
				o.close()
				o=open("output.txt",'a')
				pickle.dump(generation_moves,open("generation_moves.p",'wb'))
				if itr < iterations:
					if serial == True:
						scalc=len(generation_energies)
						if generation_complete==False:
							for i in range(0+scalc,search_gen):
								if itr == iterations:
									break			
								nstruct=len(generation)
								o.write("\n")
								atoms=generation[i][0]
								iat=len(atoms)
								print (str("structure " + str(i+1) +" of " + str(nstruct)),end="\r")

								if ctype == 'gulp':
									try:
										atoms,energy,converged=run_gulp(atoms=atoms,shel=shel,kwds=kwds,opts=gulp_opts,lib=lib,produce_steps=produce_steps,gulp_command=gulp_command,gulp_timeout=gulp_timeout)
									except:
										converged = False
										energy = 1.e20
										
								if ctype == 'vasp':
									try:
										atoms,energy,converged=run_vasp(atoms=atoms,vasp_opts=vasp_opts,kcut=kcut,produce_steps=produce_steps,dist_cutoff=dist_cutoff)
									except:
										converged = False
										energy = 1.e20	
								
								if ctype == 'qe':
									try:
										atoms,energy,converged=run_qe(atoms=atoms,qe_opts=qe_opts,kcut=kcut,produce_steps=produce_steps)
									except:
										converged=False
										energy=1.e20

								if ctype == 'mixed':		
									try:
										atoms,energy,converged=run_calculators(atoms=atoms,vasp_opts=
											vasp_opts,kcut=kcut,produce_steps=produce_steps,shel=shel,
											kwds=kwds,gulp_opts=gulp_opts,lib=lib,calcs=calcs,dist_cutoff=dist_cutoff,qe_opts=qe_opts,
											gulp_command=gulp_command,gulp_timeout=gulp_timeout)
						
									except:
										converged = False
										energy = 1.e20
		
								if len(atoms) != iat:
									converged = False
									energy = 1.e20
								if converged != True:
									if produce_steps==True:
										targets=glob.glob("atoms*.cif")
										targets.sort()
										for z in range(len(targets)):
											label=str("steps/S"+str(count)+"_"+str(z+1)+str(".cif"))
											shutil.copy(targets[z],label)

								energy = energy/len(atoms)
								energy = float(Decimal(energy).quantize(Decimal('1e-6')))
								generation[i].append(energy)
								generation[i][0]=atoms
								#print(str("I{0:=6n}".format(i+1)+"  {0:0=6.6n}	  eV/atom".format(float(energy)).rjust(25)))
								if search == 1:
									o.write(str("BH{0:=7n}".format(count)+"	 {0:=6.6n}   eV/atom".format(float(energy)).rjust(25)))
	
								if write_all_structures == True:
									try:
										write(str("structures/"+str("S-")+str('{0:0=7n}'.format(count)+".cif")),atoms)
									except:
										pass
								#if write_graph==True:
								try:
									if converged != False:
										graph['move'].append(str(generation_moves[i]))
										graph['type'].append("S")
										graph['step'].append(count+initial_gen)
										graph['energies'].append(energy)
										if search == 1:
											graph['temp'].append(T)
										graph['current_energy'].append(current_structure[3])
										graph_write=pandas.DataFrame(graph)
										graph_write.to_csv(path_or_buf="graph_output.csv",index=False)
									if count % 10 == 0:
										if write_graph==True:
												plot_graph(search=search)
								except:
									graph['move']=[str(generation_moves[i])]
									graph['type']=["S"]
									graph['step']=[count+initial_gen]
									graph['energies']=[energy]
									if search == 1:
										graph['temp']=[T]
									graph['current_energy']=[current_structure[-1]]
									graph_write=pandas.DataFrame(graph)
									graph_write.to_csv(path_or_buf="graph_output.csv",index=False)
								if count % 10 == 0:
									if write_graph==True:
										plot_graph(search=search)
								#print(graph)
									
								generation_energies.append(energy)
								search_structures[count]=[atoms,generation[i][1],generation[i][2],energy,generation_moves[i]]
								itr+=1
								count+=1
								if len(generation_energies) == search_gen:
									generation_complete=True
									print("Generation completed: ",end=' ')
									o.write("\nGeneration completed: ")
								if os.path.isfile("stop.txt"):
									print ("** Stopping caclulation **")
									### while this will now run all of the initial population, need to then log & store the data!
									# now need to write out what we've got so far...
									print ("Writing restart files...")
									o.write("\nWriting restart files...\n")
									min_structure_num=energies.index(min(energies))
									min_structure=search_structures[min_structure_num][0]
									write("lowest_energy_structure.cif",min_structure)
									pickle.dump(initial_population,open("initial_structures.p",'wb'))
									pickle.dump(energies,open("energies.p",'wb'))
									pickle.dump(search_structures,open("search_structures.p",'wb'))
									if initial_complete==1:
										pickle.dump(current_structure,open("current_structure.p",'wb'))
									pickle.dump(generation_complete,open("generation_complete.p",'wb'))
									if generation_complete==False:
										pickle.dump(generation,open("generation.p",'wb'))
										pickle.dump(generation_energies,open("generation_energies.p",'wb'))
									pickle.dump(search_structures,open("search_structures.p",'wb'))
									pickle.dump(r,open("r.p",'wb'))
									pickle.dump(ini_energies,open("ini_energies.p",'wb'))
									pickle.dump(ca,open("ca.p",'wb'))
									pickle.dump(T,open("T.p",'wb'))
									pickle.dump(search,open("search.p",'wb'))
									pickle.dump(moves,open("moves.p",'wb'))
									#if write_graph == True:
										
									graph_to_write=graph.copy()
									graph_to_write['file_name']=[]
									for x in range(len(graph_write['step'])):
										number=graph_to_write['step'][x]
										if graph_to_write['type'][x] == 'S':
											number -= initial_gen
										
										label=str( graph_to_write['type'][x] + "_" + str("{0:06d}").format(number) + ".cif")
										graph_to_write['file_name'].append(label)
									
									graph_to_write=pandas.DataFrame(graph_to_write)
									graph_to_write.to_csv(path_or_buf="graph_output.csv",index=False)
									if write_graph==True:
										plot_graph(search=search)

									#############################################################################
									
									### print out total runtime #################################################
									t2=datetime.datetime.now()
									print("\ntotal time: "+str(t2-t1)+" hours:minutes:seconds")
									o.write("\n\ntotal time: "+str(t2-t1)+" hours:minutes:seconds\n")
									o.close()	
									os.remove("stop.txt")
									sys.exit()
								else:
									min_structure_num=energies.index(min(energies))
									min_structure=search_structures[min_structure_num][0]
									write("lowest_energy_structure.cif",min_structure)
									pickle.dump(initial_population,open("initial_structures.p",'wb'))
									pickle.dump(energies,open("energies.p",'wb'))
									pickle.dump(search_structures,open("search_structures.p",'wb'))
									if initial_complete==1:
										pickle.dump(current_structure,open("current_structure.p",'wb'))
									pickle.dump(generation_complete,open("generation_complete.p",'wb'))
									if generation_complete==False:
										pickle.dump(generation,open("generation.p",'wb'))
										pickle.dump(generation_energies,open("generation_energies.p",'wb'))
									pickle.dump(search_structures,open("search_structures.p",'wb'))
									pickle.dump(r,open("r.p",'wb'))
									pickle.dump(ini_energies,open("ini_energies.p",'wb'))
									pickle.dump(ca,open("ca.p",'wb'))
									pickle.dump(T,open("T.p",'wb'))		
									pickle.dump(moves,open("moves.p",'wb'))
									pickle.dump(search,open("search.p",'wb'))
						if generation_complete == True:
							#if write_graph==True:
							#	plot_graph(search)
							# MC Accept part #############################################
							if min(generation_energies) < current_structure[3]:
							#if min(generation_energies) < min(energies):
								idx=generation_energies.index(min(generation_energies))
								current_structure=generation[idx].copy()
								new=1
								if min(generation_energies) < min(energies):
									r=0
									sa=0
								if min(generation_energies) != min(energies):
									ca=0
									if search == 1:
										moves = n_moves		
								T=T_0
								if min(generation_energies) < min(energies):
									write("lowest_energy_structure.cif",current_structure[0])
								
								write("current_structure.cif",current_structure[0])
							if min(generation_energies) >= min(energies):
								if min(generation_energies) < 0:
									new=0
									if search == 1:
										rand=random.random()
										diff=min(generation_energies)-min(energies)
										Test=math.exp(-diff/T)
										if Test>=rand:
											new=2
											idx=generation_energies.index(min(generation_energies))
											current_structure=generation[idx].copy()
											write("current_structure.cif",current_structure[0])
											
											#T=T_0
											#ca=0
								r+=len(generation)
								ca+=len(generation)
							#print(("Ca "+str(ca)),end='')	
							##############################################################
							dE=float(Decimal(min(generation_energies)-min(energies)).quantize(Decimal('1e-6')))
							print(str("E = "+str("{0: .5e}").format(min(generation_energies)) + " dE vs. global = " + str("{0: .4e}").format(dE).rjust(7)+" r: "+str(r).rjust(4)),end='')
		
							if search == 1:
								print("	 T = "+str("{0:.5f}").format(T))
							else:
								print('')
							#print(sa)
							
							o.write("E = "+str("{0: .5e}").format(min(generation_energies)) + " dE vs. global = " + str("{0: .4e}").format(dE).rjust(7)+" r: "+str(r).rjust(4))

							if search == 1:
								o.write("  T = "+str("{0:.5f}").format(T)+"\n")
							else:
								o.write("\n")		
						
							for j in range(len(generation_energies)):
								energies.append(generation_energies[j])
														
							#print(generation_energies)
							#try:
							#	start_point=max(list(search_structures.keys()))
							#except:
							#	start_point=0
							#for j in range(len(generation)):
							#	search_structures[start_point+j+1]=generation[j]		
							o.flush()
					#######################################################################
				
					#######################################################################	
					
					if search == 1:
						if ca >= melt_threshold:
							T+=(random.random()/100)
							T=float(Decimal(T).quantize(Decimal('1e-6')))
							moves=r_moves
								
				if itr >= iterations:
					break	
							
			
					#############################################################################
	
	if r >= rmax:
		print("\n*** Break conditions met stopping calculation ***\n")
		o.write("\n\n*** Break conditions met stopping calculation ***\n")
		print("Global minimum energy: "+str(min(energies))+" Structure: "+str(energies.index(min(energies))))
		o.write("\nGlobal minimum energy: "+str(min(energies))+" Structure: "+str(energies.index(min(energies))))
			
	else:
		if iterations > 0:
			print("\n*** Max number of steps reached stopping calculation ***")
			o.write("\n\n*** Max number of steps reached stopping calculation ***")
			print("Current global minimum energy: "+str(min(energies))+" Structure: "+str(energies.index(min(energies))))
			o.write("\nCurrent Global minimum energy: "+str(min(energies))+" Structure: "+str(energies.index(min(energies))))
		
	### while this will now run all of the initial population, need to then log & store the data!
	# now need to write out what we've got so far...
	print ("Writing restart files...")
	
	try:
		min_structure_num=energies.index(min(energies))
		min_structure=search_structures[min_structure_num][0]
		write("lowest_energy_structure.cif",min_structure)
	except KeyError:
		pass
	
	if write_graph_end == True:
	
		#if write_graph == True:	
		graph_to_write=graph.copy()
		graph_to_write['file_name']=[]
		for x in range(len(graph_to_write['step'])):
			number=graph_to_write['step'][x]
			if graph_to_write['type'][x] == 'S':
				number -= initial_gen
			
			label=str( graph_to_write['type'][x] + "_" + str("{0:06d}").format(number) + ".cif")
			graph_to_write['file_name'].append(label)
		
		graph_to_write=pandas.DataFrame(graph_to_write)
		graph_to_write.to_csv(path_or_buf="graph_output.csv",index=False)
		plot_graph(search=search)
		
		#else:
		#	pass
			
			
	o.write("\nWriting restart files...\n")
	pickle.dump(initial_population,open("initial_structures.p",'wb'))
	pickle.dump(energies,open("energies.p",'wb'))
	pickle.dump(search_structures,open("search_structures.p",'wb'))
	if initial_complete==1:
		pickle.dump(current_structure,open("current_structure.p",'wb'))
	pickle.dump(generation_complete,open("generation_complete.p",'wb'))
	if generation_complete==False:
		pickle.dump(generation,open("generation.p",'wb'))
		pickle.dump(generation_energies,open("generation_energies.p",'wb'))
	pickle.dump(search_structures,open("search_structures.p",'wb'))
	pickle.dump(r,open("r.p",'wb'))
	pickle.dump(ini_energies,open("ini_energies.p",'wb'))
	pickle.dump(ca,open("ca.p",'wb'))
	pickle.dump(T,open("T.p",'wb'))
	pickle.dump(search,open("search.p",'wb'))
	pickle.dump(moves,open("moves.p",'wb'))
	
	
	#############################################################################
	
	### print out total runtime #################################################
	t2=datetime.datetime.now()
	print("\ntotal time: "+str(t2-t1)+" hours:minutes:seconds")
	o.write("\n\ntotal time: "+str(t2-t1)+" hours:minutes:seconds\n")
	o.close()		
