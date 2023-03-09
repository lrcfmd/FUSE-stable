from ase import *
import pickle
from ase.visualize import *
import os
import sys
import numpy

#script to read in all of the structures from a fuse run and write them as cif files
#ideally want options to: write all: individually & as a multi-cif, or write specific
#/ranges/list of structures

def write_structures():
	try:
		ini_structures=pickle.load(open("initial_structures.p",'rb'))
		search_structures=pickle.load(open("search_structures.p",'rb'))
		#ask the user what to do
		#print("what format would you like the structures in? <cif>(one multi-cif file)")
		#format=str(input())
		possible=['cif']
		format='cif'
		if format in possible:
			if not os.path.isdir("structures"):
				os.mkdir("structures")
		
		if format == 'cif':
			keys=list(ini_structures.keys())
			for i in range(len(keys)):
				write(str("structures/I_"+str("{0:06d}").format(keys[i])+".cif"),ini_structures[keys[i]][0])
			keys=list(search_structures.keys())
			for i in range(len(keys)):
				write(str("structures/S_"+str("{0:06d}").format(keys[i])+".cif"),search_structures[keys[i]][0])
	except:
		pass			
	
write_structures()
