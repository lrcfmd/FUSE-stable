from ase import *
from ase.io import *
from ase.calculators.vasp import Vasp
import os
import math
import sys

################################################################################################
import shlex
import re
import numpy
from numpy import arccos, pi, dot
from numpy.linalg import norm

def cellpar(atoms):
	cell = atoms.cell
	a = norm(cell[0])
	b = norm(cell[1])
	c = norm(cell[2])
	alpha = arccos(dot(cell[1], cell[2])/(b*c))*180./pi
	beta  = arccos(dot(cell[0], cell[2])/(a*c))*180./pi
	gamma = arccos(dot(cell[0], cell[1])/(a*b))*180./pi

	cell = []
	cell=[a,b,c,alpha,beta,gamma]
	return cell
################################################################################################

def run_vasp(atoms='',vasp_opts='',kcut='',produce_steps='',dist_cutoff = 1.0):
	
	short_contact = False
	converged = False
	energy = 1.e20
	temp_atoms=atoms.repeat([2,2,2])
	temp1=temp_atoms.get_all_distances()
	temp2=[]
	for i in range(len(temp1)):
		for j in range(len(temp1[i])):
			if temp1[i][j] != 0:
				temp2.append(temp1[i][j])
	distances=min(temp2)
	if distances <= dist_cutoff:
		short_contact = True

	#print("\n",short_contact,"\n")
	if short_contact == False:
		new_atoms=atoms.copy()
		for i in range(len(list(vasp_opts.keys()))):
		#	print("\n\nhello!!\n\n")
   	
			if len(kcut) >= 2:
				temp_kcut=kcut[i]
			else:
				temp_kcut=kcut
			cell=new_atoms.get_cell_lengths_and_angles()
			kp=[]
			kp=(int(math.ceil(temp_kcut/cell[0])),int(math.ceil(temp_kcut/cell[1])),int(math.ceil(temp_kcut/cell[2])))
			if os.path.isfile("WAVECAR"):
				os.remove("WAVECAR")
			calc=vasp_opts[str(list(vasp_opts.keys())[i])]
			calc.set(kpts=kp)
			new_atoms.set_calculator(calc)
			energy=new_atoms.get_potential_energy()
			if produce_steps==True:
				label=str("atoms"+str(i+1)+".cif")
				write(label,atoms)
   	
		converged=new_atoms.get_calculator().converged
		
		temp_atoms=new_atoms.repeat([2,2,2])
		temp1=temp_atoms.get_all_distances()
		temp2=[]
		for i in range(len(temp1)):
			for j in range(len(temp1[i])):
				if temp1[i][j] != 0:
					temp2.append(temp1[i][j])
		if min(temp2) <= dist_cutoff :
			converged=False
   	
		
	return new_atoms,energy,converged
	

#vasp_opts={
#'1':Vasp(prec='Fast',ediff=1.E-4,nelm=50,encut=450,ibrion=2,isif=3,nsw=150,ediffg=-0.05,nwrite=1,setups={'Sr':'_sv','Ti':'_pv'},ncore=4,nelmin=4,potim=0.1,isym=0,algo='Fast',gamma=True),
#}

#atoms,energy,converged=run_vasp(atoms=read("POSCAR"),vasp_opts=vasp_opts,kcut=[20])

