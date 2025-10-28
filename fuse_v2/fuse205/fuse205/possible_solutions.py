import numpy
import pickle
import os
import sys
import platform
import glob

def return_factors(x):
   # return factors of x
   factors=[]
   for i in range(1, x + 1):
   	if x % i == 0:
   		factors.append(i)
   return factors

def cube_function(x):
	return (x**2) * (2*x)
	
def tetragonal_function(x,z):
	return (x**2) * z

def orthorhombic_function(x,w,z):
	return x*w*z


def possible_solutions(max_ax,restart):

	if restart == False:
		if platform.system == 'Windows':
			temp=glob.glob("*.p")
			for i in range(len(temp)):
				os.remove(temp[i])
		if platform.system == 'Linux':
			os.system("rm *.p")		
	
	
	#is cubic possible? #########################################################
	# cubic is defined as integer solutions to (n**2 * 2*n)-y = 0
	# below, is the total number of sub-modules required to make cubic lattices, for n = 1 - 50
	try:
		cubic_solutions=pickle.load(open("cubes.p",'rb'))
	except:
		cubic_solutions={}
		for i in range(1,max_ax+1):
			cubic_solutions[cube_function(i)]=[i,i,2*i]	
		
		pickle.dump(cubic_solutions,open("cubes.p",'wb'))	
	#############################################################################
	
	
	# possible tetragonal solutions? ############################################
	# below is a dictionary of all possible x and z pairs for tetragonal cells for y sub-modules, for x & z upto 50.
	# some numbers of sub modules have more than one solution, which is why each list is stored as a pair
	try:
		tetragonal_solutions=pickle.load(open("tetragonal.p",'rb'))
	except:
		tetragonal_solutions={}
		for i in range(1,max_ax+1):
			for j in range(2,max_ax+1,2):
				y=tetragonal_function(i,j)
				try:
					tetragonal_solutions[y].append([i,j])
				except:
					tetragonal_solutions[y]=[[i,j]]
						
		pickle.dump(tetragonal_solutions,open("tetragonal.p",'wb'))	
					
	#############################################################################
	
	# possible hexagonal solutions ##############################################
	# produced from the same as above, but without the z = even restriction 
	# some numbers of sub modules have more than one solution, which is why each list is stored as a pair
	try:
		hexagonal_solutions=pickle.load(open("hexagonal.p",'rb'))
	except:
		hexagonal_solutions={}
		for i in range(1,max_ax+1):
			for j in range(2,max_ax+1):
				y=tetragonal_function(i,j)
				try:
					hexagonal_solutions[y].append([i,j])
				except:
					hexagonal_solutions[y]=[[i,j]]
						
		pickle.dump(hexagonal_solutions,open("hexagonal.p",'wb'))	
					
	#############################################################################
	
	# possible orthorhombic solutions ###########################################
	#below are the possible orthorhombic solutions, stored as above, with the restriction
	# that z = even upto 50 sub-mods in every direction have to generate them on the
	# fly as having the full list in the file kills jedit!
	try:
		orthorhombic_solutions=pickle.load(open("orthorhombic.p",'rb'))
	except:
		orthorhombic_solutions={}
		for i in range(1,max_ax+1):
			for j in range(1,max_ax+1):
				for k in range(2,max_ax+1,2):
					y=orthorhombic_function(i,j,k)
					try:
						orthorhombic_solutions[y].append([i,j,k])
					except:
						orthorhombic_solutions[y]=[[i,j,k]]
		pickle.dump(orthorhombic_solutions,open("orthorhombic.p",'wb'))	

	#print(orthorhombic_solutions)
	
	#############################################################################
	
	# possible monoclinic & triclinic dimensions, as orthorhombic but without the
	# z = evan restriction
	try:
		monoclinic_solutions=pickle.load("monoclinic.p",'rb')
	except:
		monoclinic_solutions={}
		for i in range(1,max_ax+1):
			for j in range(1,max_ax+1):
				for k in range(1,max_ax+1):
					y=orthorhombic_function(i,j,k)
					try:
						monoclinic_solutions[y].append([i,j,k])
					except:
						monoclinic_solutions[y]=[[i,j,k]]
		pickle.dump(monoclinic_solutions,open("monoclinic.p",'wb'))	
					
	#print(orthorhombic_solutions)
	
	#############################################################################
	
	return cubic_solutions,tetragonal_solutions,hexagonal_solutions,orthorhombic_solutions,monoclinic_solutions
