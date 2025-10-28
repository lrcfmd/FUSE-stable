#from sevenn.sevennet_calculator import SevenNetCalculator
import ase
from orb_models.forcefield import pretrained
from orb_models.forcefield.calculator import ORBCalculator
import os
import glob
from ase.io import *

import sys
import time
import pandas
import statistics
import pickle

#from ase.io import read,write
from ase.visualize import *
from ase.filters import ExpCellFilter
from ase.filters import FrechetCellFilter

from ase.optimize import QuasiNewton
from ase.optimize import BFGS
from ase.optimize import MDMin
from ase.optimize import FIRE
from ase.md.langevin import Langevin

from ase.io import Trajectory

#check=os.environ['orb']

#calc=SevenNetCalculator(check)

device="cpu" # or device="cuda"
orbff = pretrained.orb_v1(device=device) # or choose another model using ORB_PRETRAINED_MODELS[model_name]()
calc = ORBCalculator(orbff)
#atoms = bulk('Cu', 'fcc', a=3.58, cubic=True)
#
#atoms.set_calculator(calc)
#atoms.get_potential_energy()

def run_orb(atoms):
	ftol=0.05
	
	try:	
		lattice,positions,numbers=spglib.standardize_cell(atoms,symprec=1.e-5)
		temp2=Atoms(numbers=numbers,pbc=True)
		temp2.cell=lattice
		temp2.set_scaled_positions(positions)
		atoms=temp2.copy()
		#f=open("spglib.log",'w')
		#print("\n\nI'm using SPGLIB! ",str(i),"\n\n")
		#f.close()
    
	except:
		#f=open("spglib.log",'w')
		#print("\n\nI failed at using SPGLIB!! ",str(i),"\n\n")
		#f.close()
		pass

	atoms.calc=calc
	ecf= FrechetCellFilter(atoms)
	qn=FIRE(ecf)
	traj=Trajectory("relax.traj","w",atoms)
	qn.attach(traj)
	qn.run(fmax=0.1,steps=500)
	
	forces=atoms.get_forces()
	if abs(forces.max()) >= 0.5:
		converged=False
		
		return atoms,energy,converged

	try:	
		lattice,positions,numbers=spglib.standardize_cell(atoms,symprec=1.e-5)
		temp2=Atoms(numbers=numbers,pbc=True)
		temp2.cell=lattice
		temp2.set_scaled_positions(positions)
		atoms=temp2.copy()
		#f=open("spglib.log",'w')
		#print("\n\nI'm using SPGLIB! ",str(i),"\n\n")
		#f.close()
    
	except:
		#f=open("spglib.log",'w')
		#print("\n\nI failed at using SPGLIB!! ",str(i),"\n\n")
		#f.close()
		pass

#	#intermediate step, run short MD run
#	T=1000 #T in kelvin
#	total_time=5000 #time in fs
#	step_size = 5 #step size in fs
#	steps=total_time/step_size
#	dyn = Langevin(atoms, step_size * units.fs, T * units.kB, 0.002)

#	def printenergy(a=atoms):  # store a reference to atoms in the definition.
#	    """Function to print the potential, kinetic and total energy."""
#	    epot = a.get_potential_energy() / len(a)
#	    ekin = a.get_kinetic_energy() / len(a)
#	    print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
#	          'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))
	
	
#	dyn.attach(printenergy, interval=50)
	
#	# We also want to save the positions of all atoms after every 100th time step.
#	traj = Trajectory('moldyn3.traj', 'w', atoms)
#	dyn.attach(traj.write, interval=50)
	
	# Now run the dynamics
#	printenergy()
#	dyn.run(steps)
	
	#part 3, relax the full structure
	ecf= FrechetCellFilter(atoms)
	qn=BFGS(ecf)
    #qn=BFGS(ecf)
	traj=Trajectory("relax.traj","w",atoms)
	qn.attach(traj)
	qn.run(fmax=ftol,steps=500)
	#write(j+"_after.cif",atoms)
	
	energy=atoms.get_potential_energy()
	
	forces=atoms.get_forces()
	if abs(forces.max()) <= ftol*1.2:
		converged=True
	else:
		converged=False
	#print("\n\n****",atoms,energy,converged,"\n\n")
	
	return atoms,energy,converged
