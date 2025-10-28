from ase import *
from ase.io import *
from ase.calculators.vasp import Vasp
from fuse205.gulp import *
from ase.calculators.gulp import GULP
from ase.calculators.castep import Castep
import platform
import re
import os
import math
import sys 
from fuse205.run_castep import *
from fuse205.run_gulp import *
from fuse205.run_vasp import *
from fuse205.run_qe import *

import glob
import spglib

################################################################################################
import shlex
import re
import numpy
from numpy import arccos, pi, dot
from numpy.linalg import norm

#def cellpar(atoms):
#	cell = atoms.cell
#	a = norm(cell[0])
#	b = norm(cell[1])
#	c = norm(cell[2])
#	alpha = arccos(dot(cell[1], cell[2])/(b*c))*180./pi
#	beta  = arccos(dot(cell[0], cell[2])/(a*c))*180./pi
#	gamma = arccos(dot(cell[0], cell[1])/(a*b))*180./pi
#
#	cell = []
#	cell=[a,b,c,alpha,beta,gamma]
#	return cell
################################################################################################
# making the assumption that we don't want the energy from gulp when using this, trying setting the 
# energy to zero after gulp has run
def run_calculators(atoms='',vasp_opts='',kcut='',produce_steps='',rel='',
	shel=None,kwds='',gulp_opts='',lib='',calcs='',dist_cutoff='',qe_opts='',
	gulp_command='gulp < gulp.gin > gulp.got',gulp_timeout='',castep_opts='',
	n_opts='',relaxer_opts='',opt_class='',opt_device='',mode='relax',
	use_spglib='',repose_seed='',repose_command='',repose_devmax=''):

	converged = None
	energy=0
	for x in range(len(calcs)):
		short_contact = False
		
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
		
		if short_contact == False:
			if energy < 50.:
				if calcs[x] == 'gulp':
					atoms, energy, converged = run_gulp(atoms=atoms,shel=shel,kwds=kwds,opts=gulp_opts,lib=lib,produce_steps=produce_steps,gulp_command=gulp_command,gulp_timeout=gulp_timeout,use_spglib=use_spglib)
					energy = 0.
					
				elif calcs[x] == 'vasp':
					atoms,energy,converged=run_vasp(atoms=atoms,vasp_opts=vasp_opts,kcut=kcut,produce_steps=produce_steps,dist_cutoff=dist_cutoff,use_spglib=use_spglib)
		
				elif calcs[x] == 'castep':
					atoms,energy,converged=run_castep(atoms=atoms,castep_opts=castep_opts,dist_cutoff=dist_cutoff,use_spglib=use_spglib)

				elif calcs[x] == 'qe':
					atoms,energy,converged=run_qe(atoms=atoms,qe_opts=qe_opts,kcut=kcut,produce_steps=produce_steps,use_spglib=use_spglib)
					
				elif calcs[x] == 'chgnet':
					from fuse205.run_chgnet import run_chgnet
					atoms,energy,converged = run_chgnet(atoms,n_opts=n_opts,rel=rel,relaxer_opts=relaxer_opts,opt_class=opt_class,opt_device=opt_device,mode=mode,use_spglib=use_spglib)

				elif calcs[x] == 'orb':
					from fuse205.run_orb import run_orb
					atoms,energy,converged = run_orb(atoms)

				elif calcs[x] == 'repose':
					try:
						write(repose_seed+".cell",atoms)	
						os.system(repose_command+" "+repose_seed+" > repose.out")
						try:
							output=open("repose.out",'r').readlines()
						except:
							output=None
						energy=None
						dev=0
						deviation=None
						if output!= None:
							for z in output:
								if "Enthalpy:" in z:
										energy=z
								if "Deviation:" in z:
										deviation=z
    						
							energy = float(energy.split(":")[-2])
							if deviation != None:
								dev = float(deviation.split(":")[-1])
							atoms=read(repose_seed+"-out.cell")
							converged=True
      		
					except:
						if not calcs[x] == 'repose':
							converged=False
							energy=1.e20
						else:
							converged=False
							energy=0.

					if calcs[x] == 'repose':	
						if dev > repose_devmax:
							converged=False
							energy=0.

	try:
		if glob.glob("gulptmp*") != []:
			if platform.system() == 'Windows':
				files=glob.glob("gulptmp*")
				for z in range(len(files)):
					os.remove(files[z])
			if platform.system() == 'Linux':
				os.system("rm gulptmp*")
	except:
		pass

	
	return atoms,energy,converged
