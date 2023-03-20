import sys
import pickle
from random import choice

def return_factors(x):
   # return factors of x
   factors=[]
   for i in range(1, x + 1):
   	if x % i == 0:
   		factors.append(i)
   return factors

### function to generate a random set of assembly instructions for a crystal structure
def create_random_instructions(string,ap,cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,monoclinic_solutions,instructions):
	
	#print(instructions)
	if string == None:
		instructions = None
		return instructions
	
	if len(instructions) == 2:
		nsub = int(len(string)/4)
		lattice=instructions[1]
		# chose a set of cell dimensions ############################################
		if lattice == 0: #cubic
			cell=cubic_solutions[nsub]
			full=[cell[0],cell[1],cell[2]]
			angles=[90,90,90]
			
		if lattice == 1: #tetragonal
			cell=choice(tetragonal_solutions[nsub])
			full=[cell[0],cell[0],cell[1]]
			angles=[90,90,90]
			
		if lattice == 2: #hexagonal
			cell=choice(hexagonal_solutions[nsub])
			full=[cell[0],cell[0],cell[1]]
			angles=[90,90,120]
			
		if lattice == 3: #orthorhombic
			cell=choice(orthorhombic_solutions[nsub])
			full=[cell[0],cell[1],cell[2]]
			angles=[90,90,90]
		
		if lattice == 4: #monoclinic
			cell=choice(monoclinic_solutions[nsub])
			full=[cell[0],cell[1],cell[2]]
			angles=[90,90,120]
		
		if lattice == 5: #triclinic
			cell=choice(monoclinic_solutions[nsub])
			full=[cell[0],cell[1],cell[2]]
			if full[2] % 2 == 0:
				angles=[choice([60,75,90,105,120,150]),choice([60,75,90,105,120,150]),choice([90,120])]
			else:
				angles=[choice([60,75,90,105,120,150]),choice([60,75,90,105,120,150]),120]
				
		# add in the cell dimensions into the instructions ##########################
		for i in range(len(full)):
			instructions.append(full[i])
		for i in range(len(angles)):
			instructions.append(angles[i])
		#############################################################################
		#print(instructions)
		
		#sys.exit()
		
	if len(instructions) != 2:
		instructions=[] 
		# 0 = ap parameter
		instructions.append(ap)
		
		# number of sub modules
		nsub = int(len(string)/4)
		# need to work out what shapes we can form this into 
		possible=[0,0,0,0,0,0] 
		# string of six flags, one for cubic, tetragonal, hexagonal, orthorhombic, monoclinic and triclinic lattices
		#	0 indicates not possible, 1 indicates possible
		# work out factors of nsub to work out which lattices are possibe
		temp=return_factors(nsub)
		
		#is cubic possible? #########################################################
		# cubic is defined as integer solutions to (n**2 * 2*n)-y = 0
		if nsub in cubic_solutions:
			possible[0] = 1
			cube_solution=cubic_solutions[nsub][0]
		#############################################################################
		
		
		# possible tetragonal solutions? ############################################
					
		if nsub in list(tetragonal_solutions.keys()):
			possible[1]=1
			tetragonal_solution=tetragonal_solutions[nsub]
		
		#############################################################################
		
		# possible hexagonal solutions ##############################################
		if nsub in list(hexagonal_solutions.keys()):
			possible[2]=1
			hexagonal_solution=hexagonal_solutions[nsub]
		
		
		#############################################################################
		
		# possible orthorhombic solutions ###########################################
		if nsub in list(orthorhombic_solutions.keys()):
			possible[3]=1
			orthorhombic_solution=orthorhombic_solutions[nsub]
		
		#############################################################################
		
		# possible monoclinic & triclinic dimensions, as orthorhombic but without the
		if nsub in list(monoclinic_solutions.keys()):
			possible[4]=1
			possible[5]=1
			monoclinic_solution=monoclinic_solutions[nsub]
		
		#############################################################################
		
		### check to see that something is possible #################################
		test=possible.count(1) == 0
		if test == True:
			accept=0
			instructions = None
			return instructions
		
		#############################################################################
		
		# from possible lattices, chose one! ########################################
		# string of six flags, one for cubic, tetragonal, hexagonal, orthorhombic, monoclinic and triclinic lattices
		num_poss=[i for i, x in enumerate(possible) if x == 1]
		#print(num_poss)
		lattice=choice(num_poss) 
		# 0 = cubic
		# 1 = tetragonal
		# 2 = hexagonal
		# 3 = orthorhombic
		# 4 = monoclinic
		# 5 = triclinic
		# append lattice choice as element 1 in instructions
		instructions.append(lattice)
		#############################################################################
		
		# chose a set of cell dimensions ############################################
		if lattice == 0: #cubic
			cell=cube_solution
			full=[cell,cell,2*cell]
			angles=[90,90,90]
			
		if lattice == 1: #tetragonal
			cell=choice(tetragonal_solution)
			full=[cell[0],cell[0],cell[1]]
			angles=[90,90,90]
			
		if lattice == 2: #hexagonal
			cell=choice(hexagonal_solution)
			full=[cell[0],cell[0],cell[1]]
			angles=[90,90,120]
			
		if lattice == 3: #orthorhombic
			cell=choice(orthorhombic_solution)
			full=[cell[0],cell[1],cell[2]]
			angles=[90,90,90]
		
		if lattice == 4: #monoclinic
			cell=choice(monoclinic_solution)
			full=[cell[0],cell[1],cell[2]]
			angles=[90,90,120]
		
		if lattice == 5: #triclinic
			cell=choice(monoclinic_solution)
			full=[cell[0],cell[1],cell[2]]
			if full[2] % 2 == 0:
				angles=[choice([60,75,90,105,120,150]),choice([60,75,90,105,120,150]),choice([90,120])]
			else:
				angles=[choice([60,75,90,105,120,150]),choice([60,75,90,105,120,150]),120]
				
			
		# add in the cell dimensions into the instructions ##########################
		for i in range(len(full)):
			instructions.append(full[i])
		for i in range(len(angles)):
			instructions.append(angles[i])
				
		#############################################################################
		
	# need to generate translation instructions #################################
	if full[2] != 1:
		if lattice == 0 or lattice == 1 or lattice == 3: # gamma = 90 lattices
			trans_pattern=[]
			for i in range(full[2]):
				if (i+1) % 2 == 0:
					trans_pattern.append(1)
				else:
					trans_pattern.append(0)
   	
		elif lattice == 2 or lattice == 4: # gamma = 120 lattices
			trans_pattern=[]
			for i in range(full[2]):
				if i == 0:
					trans_pattern.append(0)
				elif i != full[2]-1:
					c=choice([0,1,2])
					while c == trans_pattern[-1]:
						c=choice([0,1,2])
					trans_pattern.append(c)
				if i == (full[2]-1):
					temp=[0,1,2]
					try:
						temp.remove(trans_pattern[0])
					except:
						pass
					try:
						temp.remove(trans_pattern[-1])
					except:
						pass
					trans_pattern.append(choice(temp))
   
		else: # triclinic lattice 
			if angles[-1] == 90:
				trans_pattern=[]
				for i in range(full[2]):
					if (i+1) % 2 == 0:
						trans_pattern.append(1)
					else:
						trans_pattern.append(0)
						
			if angles[-1] == 120:
				trans_pattern=[]
				for i in range(full[2]):
					if i == 0:
						trans_pattern.append(0)
					elif i != full[2]-1:
						c=choice([0,1,2])
						while c == trans_pattern[-1]:
							c=choice([0,1,2])
						trans_pattern.append(c)
					if i == (full[2]-1):
						temp=[0,1,2]
						try:
							temp.remove(trans_pattern[0])
						except:
							pass
						try:
							temp.remove(trans_pattern[-1])
						except:
							pass
						trans_pattern.append(choice(temp))
	
	if full[2] == 1:
		trans_pattern=[0]
   
   
	for i in range(len(trans_pattern)):
		instructions.append(trans_pattern[i])
		
	#############################################################################
	#print(instructions)
	return instructions
