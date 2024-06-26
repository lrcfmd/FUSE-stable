import sys
from ase import *
from ase.visualize import *
from random import choice
from fuse106.assemble_structure import *
from fuse106.create_random_string import *
from fuse106.create_random_instructions import *
from fuse106.generate_random_structure import *
from random import shuffle

#string=create_random_string(atoms_per_fu,fu,vac_ratio=vac_ratio,max_fus=max_fus,system_type=system_type,composition=composition)
### ^^^ this gives the random string required to make a structure, now 
#generate random assembly instructions 
#instructions=create_random_instructions(string,ap,cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,monoclinic_solutions)


def make_basin_move(current_structure,atoms_per_fu,fu,vac_ratio,max_fus,system_type,composition,ap,cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,monoclinic_solutions,
	density_cutoff,ideal_density,moves='',max_atoms='',initial_population='',search_structures='',
	check_bonds='',btol='',bondtable='',
	check_distances='',dist_cutoff='',imax_atoms=''):
	#view(current_structure[0])
	#print(current_structure)
	complete = 0
	# create a copy of the input structure #####################################
	origial_string=current_structure[1]
	original_instructions=current_structure[2]
	original_atoms,original_instructions=assemble_structure(origial_string,original_instructions)
	# make list of moves to choose from #########################################
	keys=list(moves.keys())
	m=[]
	for i in range(len(keys)):
		for j in range(moves[keys[i]]):
			m.append(keys[i])
	#############################################################################
	
	# need to start by splitting out the module string & instructions ###########
	string=current_structure[1]
	instructions=current_structure[2]
	translations=instructions[8:]
	#############################################################################
	attempts = 0
	max_attempts=100
	while complete == 0:
		if attempts == 0:
			move=choice(m)
			
		if attempts >= max_attempts:
			#print("reset move")
			move=choice(m)
			attempts=0
			
		attempts+=1
		#print(move)
		
		string=origial_string.copy()
		atoms,instructions=assemble_structure(string,original_instructions)
		### first basic move, swap two elements in the string object, original type 1,
		### 1/3 : 1/3 : 1/3 chances to swap 2 elements, swap multiple elements, totally randomise - original choice
		### now trying to heavily bias this towrads swap 2
		if move == 1:
			if system_type=='neutral':
				#view(atoms)
				seti=[0,0,0,0,0,0,0,0,1,1,1,2]
				setting=choice(seti)
				setting = 0 #***************
				if setting==0: #swap 2
					a=choice(list(range(len(string))))
					b=choice(list(range(len(string))))
						
					while string[a] == string[b]:
						b=choice(list(range(len(string))))
						
					a1=string[a]
					b1=string[b]
					
					string[a]=b1
					string[b]=a1
					try:
						atoms,instructions=assemble_structure(string,instructions)
					#view(atoms)
					except: 
						continue
					accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
					if accept == 1:
						new_structure=[atoms,string,instructions]
						complete = 1
				
				if setting == 1: #swap multiple
					to_swap=choice(list(range(len(string))))
					#print(to_swap)
					indicies=[]
					elements=[]
					for i in range(to_swap):
						target=choice(list(range(len(string))))
						while target in indicies:
							target=choice(list(range(len(string))))
						indicies.append(target)
						
					for i in range(len(indicies)):
						elements.append(string[indicies[i]])
					
					#print(elements)
					shuffle(elements)
					#print(elements)
					for i in range(len(indicies)):
						string[indicies[i]]=elements[i]
					try:
						atoms,instructions=assemble_structure(string,instructions)
					except:
						continue
					#view(atoms)
					#sys.exit()
					
					accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
					if accept == 1:
						new_structure=[atoms,string,instructions]
						complete=1
				
				if setting == 2: # totally randomise
					shuffle(string)
					try:
						atoms,instructions=assemble_structure(string,instructions)
					except:
						continue
					
					accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
					if accept == 1:
						new_structure=[atoms,string,instructions]
						complete=1
			
			if system_type == 'ionic':
				#view(atoms)
				seti=[0,0,0,0,0,0,0,0,1,1,1,2]
				setting1=choice(seti) # 0 swap 2 atoms, 1 swap random number, 2 swap all
				setting1=0 # *********************
				seti=[0,0,0,0,1,1,1,1,2] # again now want to try to bias this towards keeping groups together
				setting2=choice(seti) # 0 swap cations, 1 swap anions, 2 swap all
				seti=[0,1] # 0 ions only, 1 include vacencies
				setting3=choice(seti)
				#print(string) #### remove when completed function
				#setting = 1
				cations=[]
				anions=[]
				keys=list(composition.keys())
				for i in range(len(keys)):
					if composition[keys[i]][-1] > 0:
						a=Atoms(keys[i])
						if not a.get_atomic_numbers()[0] in cations:
							cations.append(a.get_atomic_numbers()[0])
						if setting3 ==1:
							cations.append(120)
							
					if composition[keys[i]][-1] < 0:
						a=Atoms(keys[i])
						if not a.get_atomic_numbers()[0] in anions:
							anions.append(a.get_atomic_numbers()[0])
						if setting3 == 1:
							anions.append(120)		
				
				if setting1==0: #swap 2
					if setting2 == 0:
						cation_pairs=[]
						for i in range(len(string)):
							if string[i] in cations:
								cation_pairs.append([string[i],i])
						
						if len(cations) >= 2:
							a=choice(cation_pairs)
							b=choice(cation_pairs)
							while a == b:
								b=choice(cation_pairs)
							
							string[b[1]]=a[0]
							string[a[1]]=b[0]
							try:
								atoms,instructions=assemble_structure(string,instructions)
							except:
								continue
								
							accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
	
							if accept == 1:
								new_structure=[atoms,string,instructions]
								complete=1

					if setting2 == 1:
						anion_pairs=[]
						for i in range(len(string)):
							if string[i] in anions:
								anion_pairs.append([string[i],i])
						
						if len(anions) >= 2:
							a=choice(anion_pairs)
							b=choice(anion_pairs)
							while a == b:
								b=choice(anion_pairs)
							
							string[b[1]]=a[0]
							string[a[1]]=b[0]
							try:
								atoms,instructions=assemble_structure(string,instructions)
							except:
								continue
								
							accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
							if accept == 1:
	
								new_structure=[atoms,string,instructions]
								complete=1

						
					if setting2 == 2:
						a=choice(list(range(len(string))))
						b=choice(list(range(len(string))))
							
						while string[a] == string[b]:
							b=choice(list(range(len(string))))
							
						a1=string[a]
						b1=string[b]
						
						string[a]=b1
						string[b]=a1
						try:
							atoms,instructions=assemble_structure(string,instructions)
						#view(atoms)
						except: 
							continue
						accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
						if accept == 1:
							new_structure=[atoms,string,instructions]
							complete = 1
				
				if setting1 == 1: #swap multiple
					
					if setting2 == 0:
						cation_pairs=[]
						for i in range(len(string)):
							if string[i] in cations:
								cation_pairs.append([string[i],i])
						
						if len(cations) >= 2:
							to_swap=choice(list(range(1,len(cation_pairs))))
							indicies=[]
							elements=[]
							for i in range(to_swap):
								target=choice(cation_pairs)
								while target[1] in indicies:
									target=choice(cation_pairs)
									
								indicies.append(target[1])
								elements.append(target[0])
							
							shuffle(elements)
							for i in range(len(indicies)):
								string[indicies[i]]=elements[i]
							try:
								atoms,instructions=assemble_structure(string,instructions)
							except:
								continue
								
							accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
							if accept == 1:
								new_structure=[atoms,string,instructions]
								complete=1

					if setting2 == 1:
						anion_pairs=[]
						for i in range(len(string)):
							if string[i] in anions:
								anion_pairs.append([string[i],i])
						
						if len(anions) >= 2:
							to_swap=choice(list(range(1,len(anion_pairs))))
							indicies=[]
							elements=[]
							for i in range(to_swap):
								target=choice(anion_pairs)
								while target[1] in indicies:
									target=choice(anion_pairs)
									
								indicies.append(target[1])
								elements.append(target[0])
							
							shuffle(elements)
							for i in range(len(indicies)):
								string[indicies[i]]=elements[i]
							try:
								atoms,instructions=assemble_structure(string,instructions)
							except:
								continue
							accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
							if accept == 1:
								new_structure=[atoms,string,instructions]
								complete=1

					
					if setting2 == 2:
						to_swap=choice(list(range(2,len(string))))
						#print(to_swap)
						indicies=[]
						elements=[]
						for i in range(to_swap):
							target=choice(list(range(len(string))))
							while target in indicies:
								target=choice(list(range(len(string))))
							indicies.append(target)
							
						for i in range(len(indicies)):
							elements.append(string[indicies[i]])
						#print(elements)
						shuffle(elements)
						#print(elements)
						for i in range(len(indicies)):
							string[indicies[i]]=elements[i]
						try:
							atoms,instructions=assemble_structure(string,instructions)
						except:
							continue
						#view(atoms)
						#sys.exit()
						accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
						if accept == 1:
							new_structure=[atoms,string,instructions]
							complete=1
				
				if setting1 == 2: # totally randomise
					
					if setting2 == 0:
						cation_pairs=[]
						for i in range(len(string)):
							if string[i] in cations:
								cation_pairs.append([string[i],i])
								
						if len(cations) >=2 :
							indicies=[]
							elements=[]
							for i in range(len(cation_pairs)):
								indicies.append(cation_pairs[i][1])
								elements.append(cation_pairs[i][0])
							shuffle(elements)
							for i in range(len(indicies)):
								string[indicies[i]]=elements[i]
							try:
								atoms,instructions=assemble_structure(string,instructions)
							except:
								continue
							
							accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
							if accept == 1:
								new_structure=[atoms,string,instructions]
								complete=1

					if setting2 == 1:
						anion_pairs=[]
						for i in range(len(string)):
							if string[i] in anions:
								anion_pairs.append([string[i],i])
								
						if len(anions) >=2 :
							indicies=[]
							elements=[]
							for i in range(len(anion_pairs)):
								indicies.append(anion_pairs[i][1])
								elements.append(anion_pairs[i][0])
							shuffle(elements)
							for i in range(len(indicies)):
								string[indicies[i]]=elements[i]
							try:
								atoms,instructions=assemble_structure(string,instructions)
							except:
								continue
								
							accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
							if accept == 1:
								new_structure=[atoms,string,instructions]
								complete=1

					
					if setting2 == 2:
						shuffle(string)
						try:
							atoms,instructions=assemble_structure(string,instructions)
						except:
							continue
							
						accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
						if accept == 1:
							new_structure=[atoms,string,instructions]
							complete=1			
		##########################################################################
		
		### remove two sub-modules and swap their positions ######################
		if move == 2:
			temp_string=string.copy()
			nsub=int(len(temp_string)/4)
			sub_mods=[]
			for i in range(nsub):
				temp=temp_string[-4:]
				sub_mods.append(temp)
				del temp_string[-4:]
			if len(sub_mods) > 1:
				a1=choice(list(range(len(sub_mods))))
			else:
				continue
			a2=choice(list(range(len(sub_mods))))
			while a2 == a1:
				a2=choice(list(range(len(sub_mods))))
			a1p=sub_mods[a1]
			a2p=sub_mods[a2]
			sub_mods[a1]=a2p
			sub_mods[a2]=a1p
			new_string=[]
			for i in range(len(sub_mods)):
				for j in range(len(sub_mods[i])):
					new_string.append(sub_mods[i][j])
			string=new_string
			try:
				atoms,instructions=assemble_structure(string,instructions)
				accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
				if accept == 1:
					new_structure=[atoms,string,instructions]
					complete=1
			except:
				continue
	#############################################################################
		
		### generate a new set of random instructions without changing the string
		if move == 3:
			instructions=[]
			instructions=create_random_instructions(string,ap,cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,monoclinic_solutions,instructions)
			if instructions == None:
				continue
			#send string & instructions to the structure assembler
			try:
				atoms,instructions=assemble_structure(string,instructions)
				accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
				if accept == 1:
					new_structure=[atoms,string,instructions]
					complete=1
			except:
				continue
		##########################################################################
			
		### swap the positions of two full modules ###############################
		if move == 4:
			layers=instructions[4]
			if layers > 2:
				temp_string=string.copy()
				nsub=int(len(temp_string)/4)
				sub_mods=[]
				for i in range(nsub):
					temp=temp_string[-4:]
					sub_mods.append(temp)
					del temp_string[-4:]
				sub_mods_per_layer=int(len(sub_mods)/layers)
				modules=[]
				nmod=0
				for i in range(layers):
					size=instructions[3]*instructions[2]
					layer=[]
					for j in range(int(size)):
						layer.append(sub_mods[nmod])
						nmod+=1
					modules.append(layer)
				a1=choice(list(range(len(modules))))
				b1=choice(list(range(len(modules))))
				while b1 == a1:
						b1=choice(list(range(len(modules))))

				a2=modules[a1]
				b2=modules[b1]
				modules[a1]=b2
				modules[b1]=a2
				new_string=[]
				for i in range(len(modules)):
					for j in range(len(modules[i])):
						for k in range(len(modules[i][j])):
							new_string.append(modules[i][j][k])
				string=new_string
				try:
					atoms,instructions=assemble_structure(string,instructions)
					accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
					if accept == 1:
						new_structure=[atoms,string,instructions]
						complete=1
				except:
					continue
		##########################################################################
		
		### now need to grow the structure along a random axis ###################
		### still need to write this for the a and b directions ##################
		if move == 5:
			#work out if we are currently at half or less of max atoms
			try:
				temp_atoms,instructions=assemble_structure(string,instructions)
				n_atoms=len(temp_atoms)
				factor=n_atoms/max_atoms
				if factor <= 0.5:
					# choose axis
					axis=choice([0,1,2])
					if axis == 0: # a-direction
						temp_string=string.copy()
						temp_instructions=instructions.copy()
						temp_atoms,temp_instructions2=assemble_structure(temp_string,temp_instructions)
						# extract the sub mods from the string
						new_string=[]
						sub_mods=[]
						nsub = int(len(string)/4)
						for i in range(nsub):
							temp=temp_string[-4:]
							sub_mods.append(temp)
							del temp_string[-4:]
							
						sub_mods.reverse()
						
						#now need to group sub mods into modules expanding along the desired direction
						mods=[]
						n=0
						for k in range(temp_instructions[4]):
							for i in range(temp_instructions[2]):
								tempn=n
								for j in range(temp_instructions[3]):
									mods.append(sub_mods[tempn])
									tempn+=1
								for j in range(temp_instructions[3]):
									mods.append(sub_mods[n])
									n+=1
						#now take the mods and collapse it onto a new string
						new_string=[]
						for i in range(len(mods)):
							for j in range(len(mods[i])):
								new_string.append(mods[i][j])
						
						#create the new instructions object
						new_instructions=instructions.copy()
						new_instructions[2] *= 2
						atoms,new_instructions=assemble_structure(new_string,new_instructions)
						accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
						if accept == 1:
							new_structure=[atoms,new_string,new_instructions]
							complete=1
					
					if axis == 1: # b-direction
						temp_string=string.copy()
						temp_instructions=instructions.copy()
						temp_atoms,temp_instructions2=assemble_structure(temp_string,temp_instructions)
						# extract the sub mods from the string
						new_string=[]
						sub_mods=[]
						nsub = int(len(string)/4)
						for i in range(nsub):
							temp=temp_string[-4:]
							sub_mods.append(temp)
							del temp_string[-4:]
						
						sub_mods.reverse()
						
						#now need to group sub mods into modules expanding along the desired direction
						mods=[]
						n=0
						for k in range(temp_instructions[4]):
							for i in range(temp_instructions[3]):
								for j in range(temp_instructions[2]):
									mods.append(sub_mods[n])
									mods.append(sub_mods[n])
									n+=1
						#now take the mods and collapse it onto a new string
						new_string=[]
						for i in range(len(mods)):
							for j in range(len(mods[i])):
								new_string.append(mods[i][j])
						
						#create the new instructions object
						new_instructions=instructions.copy()
						new_instructions[3] *= 2
						atoms,new_instructions=assemble_structure(new_string,new_instructions)
						accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
						if accept == 1:
							new_structure=[atoms,new_string,new_instructions]
							complete=1
					
					if axis == 2: # c-direction
						new_string=[]
						new_translations=[]
						for i in range(2):
							for j in range(len(string)):
								new_string.append(string[j])
							for j in range(len(translations)):
								new_translations.append(translations[j])
								
						new_instructions=instructions[:8]
						for i in range(len(new_translations)):
							new_instructions.append(new_translations[i])
						string=new_string
						instructions=new_instructions
						instructions[4] *= 2
						atoms,instructions=assemble_structure(string,instructions)
						accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
						if accept == 1:
							new_structure=[atoms,string,instructions]
							complete=1
			except:
				continue		
		##########################################################################

		###random new structure with similar volume: max_fus = current number of fus
		if move == 6:
			try:
				#atoms,instructions=assemble_structure(string,instructions)
				nfu=int(len(atoms)/atoms_per_fu)
				#string,instructions=create_random_string(cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,monoclinic_solutions,atoms_per_fu,fu,vac_ratio=vac_ratio,max_fus=nfu,system_type=system_type,composition=composition,ap=ap)
				#### ^^^ this gives the random string required to make a structure, now 
				##generate random assembly instructions 
				#instructions=create_random_instructions(string,ap,cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,monoclinic_solutions,instructions)
				#if instructions == None:
				#	continue
				target_fus=choice(list(range(1,nfu+1)))
				target_atoms=target_fus*atoms_per_fu
				atoms,string,instructions,accept=generate_random_structure(target_atoms,cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,
				monoclinic_solutions,atoms_per_fu,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,
				vac_ratio=vac_ratio,max_fus=max_fus,
				target_number_atoms=imax_atoms)				
				
					
				#send string & instructions to the structure assembler
				#atoms,instructions=assemble_structure(string,instructions)
				accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
				if accept == 1:
					new_structure=[atoms,string,instructions]
					complete=1
			except:
				continue
		##########################################################################
		
		### last move: random new structure ######################################
		if move == 7:
			try:
				#string,instructions=create_random_string(cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,monoclinic_solutions,atoms_per_fu,fu,vac_ratio=vac_ratio,max_fus=max_fus,system_type=system_type,composition=composition,ap=ap)
				#### ^^^ this gives the random string required to make a structure, now 
				##generate random assembly instructions 
				#instructions=create_random_instructions(string,ap,cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,monoclinic_solutions,instructions)
				#if instructions == None:
				#	continue
				#send string & instructions to the structure assembler
				#atoms,instructions=assemble_structure(string,instructions)
				target_fus=choice(list(range(1,max_fus+1)))
				target_atoms=target_fus*atoms_per_fu
				atoms,string,instructions,accept=generate_random_structure(target_atoms,cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,
				monoclinic_solutions,atoms_per_fu,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,
				vac_ratio=vac_ratio,max_fus=max_fus,
				target_number_atoms=imax_atoms)					
				
				accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
				if accept == 1:
					new_structure=[atoms,string,instructions]
					complete=1
			except:
				continue
		##########################################################################
		
		### clone of move 1, where we force it to swap more than two atoms #################################
		
		if move == 9:
			if system_type=='neutral':
				#view(atoms)
				seti=[0,0,0,0,0,0,0,0,1,1,1,2]
				setting=choice(seti)
				setting = 1 #***************
				if setting==1: #swap 2
					a=choice(list(range(len(string))))
					b=choice(list(range(len(string))))
						
					while string[a] == string[b]:
						b=choice(list(range(len(string))))
						
					a1=string[a]
					b1=string[b]
					
					string[a]=b1
					string[b]=a1
					try:
						atoms,instructions=assemble_structure(string,instructions)
					#view(atoms)
					except: 
						continue
					accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
					if accept == 1:
						new_structure=[atoms,string,instructions]
						complete = 1
				
				if setting == 1: #swap multiple
					to_swap=choice(list(range(len(string))))
					#print(to_swap)
					indicies=[]
					elements=[]
					for i in range(to_swap):
						target=choice(list(range(len(string))))
						while target in indicies:
							target=choice(list(range(len(string))))
						indicies.append(target)
						
					for i in range(len(indicies)):
						elements.append(string[indicies[i]])
					
					#print(elements)
					shuffle(elements)
					#print(elements)
					for i in range(len(indicies)):
						string[indicies[i]]=elements[i]
					try:
						atoms,instructions=assemble_structure(string,instructions)
					except:
						continue
					#view(atoms)
					#sys.exit()
					
					accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
					if accept == 1:
						new_structure=[atoms,string,instructions]
						complete=1
				
				if setting == 2: # totally randomise
					shuffle(string)
					try:
						atoms,instructions=assemble_structure(string,instructions)
					except:
						continue
					
					accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
					if accept == 1:
						new_structure=[atoms,string,instructions]
						complete=1
			
			if system_type == 'ionic':
				#view(atoms)
				seti=[0,0,0,0,0,0,0,0,1,1,1,2]
				setting1=choice(seti) # 0 swap 2 atoms, 1 swap random number, 2 swap all
				setting1=1 # *********************
				seti=[0,0,0,0,1,1,1,1,2] # again now want to try to bias this towards keeping groups together
				setting2=choice(seti) # 0 swap cations, 1 swap anions, 2 swap all
				seti=[0,1] # 0 ions only, 1 include vacencies
				setting3=choice(seti)
				#print(string) #### remove when completed function
				#setting = 1
				cations=[]
				anions=[]
				keys=list(composition.keys())
				for i in range(len(keys)):
					if composition[keys[i]][-1] > 0:
						a=Atoms(keys[i])
						if not a.get_atomic_numbers()[0] in cations:
							cations.append(a.get_atomic_numbers()[0])
						if setting3 ==1:
							cations.append(120)
							
					if composition[keys[i]][-1] < 0:
						a=Atoms(keys[i])
						if not a.get_atomic_numbers()[0] in anions:
							anions.append(a.get_atomic_numbers()[0])
						if setting3 == 1:
							anions.append(120)		
				
				if setting1==0: #swap 2
					if setting2 == 0:
						cation_pairs=[]
						for i in range(len(string)):
							if string[i] in cations:
								cation_pairs.append([string[i],i])
						
						if len(cations) >= 2:
							a=choice(cation_pairs)
							b=choice(cation_pairs)
							while a == b:
								b=choice(cation_pairs)
							
							string[b[1]]=a[0]
							string[a[1]]=b[0]
							try:
								atoms,instructions=assemble_structure(string,instructions)
							except:
								continue
								
							accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
	
							if accept == 1:
								new_structure=[atoms,string,instructions]
								complete=1

					if setting2 == 1:
						anion_pairs=[]
						for i in range(len(string)):
							if string[i] in anions:
								anion_pairs.append([string[i],i])
						
						if len(anions) >= 2:
							a=choice(anion_pairs)
							b=choice(anion_pairs)
							while a == b:
								b=choice(anion_pairs)
							
							string[b[1]]=a[0]
							string[a[1]]=b[0]
							try:
								atoms,instructions=assemble_structure(string,instructions)
							except:
								continue
								
							accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
							if accept == 1:
	
								new_structure=[atoms,string,instructions]
								complete=1

						
					if setting2 == 2:
						a=choice(list(range(len(string))))
						b=choice(list(range(len(string))))
							
						while string[a] == string[b]:
							b=choice(list(range(len(string))))
							
						a1=string[a]
						b1=string[b]
						
						string[a]=b1
						string[b]=a1
						try:
							atoms,instructions=assemble_structure(string,instructions)
						#view(atoms)
						except: 
							continue
						accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
						if accept == 1:
							new_structure=[atoms,string,instructions]
							complete = 1
				
				if setting1 == 1: #swap multiple
					
					if setting2 == 0:
						cation_pairs=[]
						for i in range(len(string)):
							if string[i] in cations:
								cation_pairs.append([string[i],i])
						
						if len(cations) >= 2:
							to_swap=choice(list(range(1,len(cation_pairs))))
							indicies=[]
							elements=[]
							for i in range(to_swap):
								target=choice(cation_pairs)
								while target[1] in indicies:
									target=choice(cation_pairs)
									
								indicies.append(target[1])
								elements.append(target[0])
							
							shuffle(elements)
							for i in range(len(indicies)):
								string[indicies[i]]=elements[i]
							try:
								atoms,instructions=assemble_structure(string,instructions)
							except:
								continue
								
							accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
							if accept == 1:
								new_structure=[atoms,string,instructions]
								complete=1

					if setting2 == 1:
						anion_pairs=[]
						for i in range(len(string)):
							if string[i] in anions:
								anion_pairs.append([string[i],i])
						
						if len(anions) >= 2:
							to_swap=choice(list(range(1,len(anion_pairs))))
							indicies=[]
							elements=[]
							for i in range(to_swap):
								target=choice(anion_pairs)
								while target[1] in indicies:
									target=choice(anion_pairs)
									
								indicies.append(target[1])
								elements.append(target[0])
							
							shuffle(elements)
							for i in range(len(indicies)):
								string[indicies[i]]=elements[i]
							try:
								atoms,instructions=assemble_structure(string,instructions)
							except:
								continue
							accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
							if accept == 1:
								new_structure=[atoms,string,instructions]
								complete=1

					
					if setting2 == 2:
						to_swap=choice(list(range(2,len(string))))
						#print(to_swap)
						indicies=[]
						elements=[]
						for i in range(to_swap):
							target=choice(list(range(len(string))))
							while target in indicies:
								target=choice(list(range(len(string))))
							indicies.append(target)
							
						for i in range(len(indicies)):
							elements.append(string[indicies[i]])
						#print(elements)
						shuffle(elements)
						#print(elements)
						for i in range(len(indicies)):
							string[indicies[i]]=elements[i]
						try:
							atoms,instructions=assemble_structure(string,instructions)
						except:
							continue
						#view(atoms)
						#sys.exit()
						accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
						if accept == 1:
							new_structure=[atoms,string,instructions]
							complete=1
				
				if setting1 == 2: # totally randomise
					
					if setting2 == 0:
						cation_pairs=[]
						for i in range(len(string)):
							if string[i] in cations:
								cation_pairs.append([string[i],i])
								
						if len(cations) >=2 :
							indicies=[]
							elements=[]
							for i in range(len(cation_pairs)):
								indicies.append(cation_pairs[i][1])
								elements.append(cation_pairs[i][0])
							shuffle(elements)
							for i in range(len(indicies)):
								string[indicies[i]]=elements[i]
							try:
								atoms,instructions=assemble_structure(string,instructions)
							except:
								continue
							
							accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
							if accept == 1:
								new_structure=[atoms,string,instructions]
								complete=1

					if setting2 == 1:
						anion_pairs=[]
						for i in range(len(string)):
							if string[i] in anions:
								anion_pairs.append([string[i],i])
								
						if len(anions) >=2 :
							indicies=[]
							elements=[]
							for i in range(len(anion_pairs)):
								indicies.append(anion_pairs[i][1])
								elements.append(anion_pairs[i][0])
							shuffle(elements)
							for i in range(len(indicies)):
								string[indicies[i]]=elements[i]
							try:
								atoms,instructions=assemble_structure(string,instructions)
							except:
								continue
								
							accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
							if accept == 1:
								new_structure=[atoms,string,instructions]
								complete=1

					
					if setting2 == 2:
						shuffle(string)
						try:
							atoms,instructions=assemble_structure(string,instructions)
						except:
							continue
							
						accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
						if accept == 1:
							new_structure=[atoms,string,instructions]
							complete=1			
		##########################################################################
		### clone of move 1, where we force it to swap all atoms #################################

		if move == 10:
			if system_type=='neutral':
				#view(atoms)
				seti=[0,0,0,0,0,0,0,0,1,1,1,2]
				setting=choice(seti)
				setting = 2 #***************
				if setting==0: #swap 2
					a=choice(list(range(len(string))))
					b=choice(list(range(len(string))))
						
					while string[a] == string[b]:
						b=choice(list(range(len(string))))
						
					a1=string[a]
					b1=string[b]
					
					string[a]=b1
					string[b]=a1
					try:
						atoms,instructions=assemble_structure(string,instructions)
					#view(atoms)
					except: 
						continue
					accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
					if accept == 1:
						new_structure=[atoms,string,instructions]
						complete = 1
				
				if setting == 1: #swap multiple
					to_swap=choice(list(range(len(string))))
					#print(to_swap)
					indicies=[]
					elements=[]
					for i in range(to_swap):
						target=choice(list(range(len(string))))
						while target in indicies:
							target=choice(list(range(len(string))))
						indicies.append(target)
						
					for i in range(len(indicies)):
						elements.append(string[indicies[i]])
					
					#print(elements)
					shuffle(elements)
					#print(elements)
					for i in range(len(indicies)):
						string[indicies[i]]=elements[i]
					try:
						atoms,instructions=assemble_structure(string,instructions)
					except:
						continue
					#view(atoms)
					#sys.exit()
					
					accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
					if accept == 1:
						new_structure=[atoms,string,instructions]
						complete=1
				
				if setting == 2: # totally randomise
					shuffle(string)
					try:
						atoms,instructions=assemble_structure(string,instructions)
					except:
						continue
					
					accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
					if accept == 1:
						new_structure=[atoms,string,instructions]
						complete=1
			
			if system_type == 'ionic':
				#view(atoms)
				seti=[0,0,0,0,0,0,0,0,1,1,1,2]
				setting1=choice(seti) # 0 swap 2 atoms, 1 swap random number, 2 swap all
				setting1=2 # *********************
				seti=[0,0,0,0,1,1,1,1,2] # again now want to try to bias this towards keeping groups together
				setting2=choice(seti) # 0 swap cations, 1 swap anions, 2 swap all
				seti=[0,1] # 0 ions only, 1 include vacencies
				setting3=choice(seti)
				#print(string) #### remove when completed function
				#setting = 1
				cations=[]
				anions=[]
				keys=list(composition.keys())
				for i in range(len(keys)):
					if composition[keys[i]][-1] > 0:
						a=Atoms(keys[i])
						if not a.get_atomic_numbers()[0] in cations:
							cations.append(a.get_atomic_numbers()[0])
						if setting3 ==1:
							cations.append(120)
							
					if composition[keys[i]][-1] < 0:
						a=Atoms(keys[i])
						if not a.get_atomic_numbers()[0] in anions:
							anions.append(a.get_atomic_numbers()[0])
						if setting3 == 1:
							anions.append(120)		
				
				if setting1==0: #swap 2
					if setting2 == 0:
						cation_pairs=[]
						for i in range(len(string)):
							if string[i] in cations:
								cation_pairs.append([string[i],i])
						
						if len(cations) >= 2:
							a=choice(cation_pairs)
							b=choice(cation_pairs)
							while a == b:
								b=choice(cation_pairs)
							
							string[b[1]]=a[0]
							string[a[1]]=b[0]
							try:
								atoms,instructions=assemble_structure(string,instructions)
							except:
								continue
								
							accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
	
							if accept == 1:
								new_structure=[atoms,string,instructions]
								complete=1

					if setting2 == 1:
						anion_pairs=[]
						for i in range(len(string)):
							if string[i] in anions:
								anion_pairs.append([string[i],i])
						
						if len(anions) >= 2:
							a=choice(anion_pairs)
							b=choice(anion_pairs)
							while a == b:
								b=choice(anion_pairs)
							
							string[b[1]]=a[0]
							string[a[1]]=b[0]
							try:
								atoms,instructions=assemble_structure(string,instructions)
							except:
								continue
								
							accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
							if accept == 1:
	
								new_structure=[atoms,string,instructions]
								complete=1

						
					if setting2 == 2:
						a=choice(list(range(len(string))))
						b=choice(list(range(len(string))))
							
						while string[a] == string[b]:
							b=choice(list(range(len(string))))
							
						a1=string[a]
						b1=string[b]
						
						string[a]=b1
						string[b]=a1
						try:
							atoms,instructions=assemble_structure(string,instructions)
						#view(atoms)
						except: 
							continue
						accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
						if accept == 1:
							new_structure=[atoms,string,instructions]
							complete = 1
				
				if setting1 == 1: #swap multiple
					
					if setting2 == 0:
						cation_pairs=[]
						for i in range(len(string)):
							if string[i] in cations:
								cation_pairs.append([string[i],i])
						
						if len(cations) >= 2:
							to_swap=choice(list(range(1,len(cation_pairs))))
							indicies=[]
							elements=[]
							for i in range(to_swap):
								target=choice(cation_pairs)
								while target[1] in indicies:
									target=choice(cation_pairs)
									
								indicies.append(target[1])
								elements.append(target[0])
							
							shuffle(elements)
							for i in range(len(indicies)):
								string[indicies[i]]=elements[i]
							try:
								atoms,instructions=assemble_structure(string,instructions)
							except:
								continue
								
							accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
							if accept == 1:
								new_structure=[atoms,string,instructions]
								complete=1

					if setting2 == 1:
						anion_pairs=[]
						for i in range(len(string)):
							if string[i] in anions:
								anion_pairs.append([string[i],i])
						
						if len(anions) >= 2:
							to_swap=choice(list(range(1,len(anion_pairs))))
							indicies=[]
							elements=[]
							for i in range(to_swap):
								target=choice(anion_pairs)
								while target[1] in indicies:
									target=choice(anion_pairs)
									
								indicies.append(target[1])
								elements.append(target[0])
							
							shuffle(elements)
							for i in range(len(indicies)):
								string[indicies[i]]=elements[i]
							try:
								atoms,instructions=assemble_structure(string,instructions)
							except:
								continue
							accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
							if accept == 1:
								new_structure=[atoms,string,instructions]
								complete=1

					
					if setting2 == 2:
						to_swap=choice(list(range(2,len(string))))
						#print(to_swap)
						indicies=[]
						elements=[]
						for i in range(to_swap):
							target=choice(list(range(len(string))))
							while target in indicies:
								target=choice(list(range(len(string))))
							indicies.append(target)
							
						for i in range(len(indicies)):
							elements.append(string[indicies[i]])
						#print(elements)
						shuffle(elements)
						#print(elements)
						for i in range(len(indicies)):
							string[indicies[i]]=elements[i]
						try:
							atoms,instructions=assemble_structure(string,instructions)
						except:
							continue
						#view(atoms)
						#sys.exit()
						accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
						if accept == 1:
							new_structure=[atoms,string,instructions]
							complete=1
				
				if setting1 == 2: # totally randomise
					
					if setting2 == 0:
						cation_pairs=[]
						for i in range(len(string)):
							if string[i] in cations:
								cation_pairs.append([string[i],i])
								
						if len(cations) >=2 :
							indicies=[]
							elements=[]
							for i in range(len(cation_pairs)):
								indicies.append(cation_pairs[i][1])
								elements.append(cation_pairs[i][0])
							shuffle(elements)
							for i in range(len(indicies)):
								string[indicies[i]]=elements[i]
							try:
								atoms,instructions=assemble_structure(string,instructions)
							except:
								continue
							
							accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
							if accept == 1:
								new_structure=[atoms,string,instructions]
								complete=1

					if setting2 == 1:
						anion_pairs=[]
						for i in range(len(string)):
							if string[i] in anions:
								anion_pairs.append([string[i],i])
								
						if len(anions) >=2 :
							indicies=[]
							elements=[]
							for i in range(len(anion_pairs)):
								indicies.append(anion_pairs[i][1])
								elements.append(anion_pairs[i][0])
							shuffle(elements)
							for i in range(len(indicies)):
								string[indicies[i]]=elements[i]
							try:
								atoms,instructions=assemble_structure(string,instructions)
							except:
								continue
								
							accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
							if accept == 1:
								new_structure=[atoms,string,instructions]
								complete=1

					
					if setting2 == 2:
						shuffle(string)
						try:
							atoms,instructions=assemble_structure(string,instructions)
						except:
							continue
							
						accept=error_check_structure(atoms,ideal_density,density_cutoff,check_bonds,btol,system_type,fu,composition,bondtable,ap,check_distances,dist_cutoff,target_number_atoms=max_atoms)
						if accept == 1:
							new_structure=[atoms,string,instructions]
							complete=1			
		##########################################################################


		

	new_structure.append(move)
		
	return new_structure
	

