from ase.io import read,write

import spglib

from ase import Atoms

from chgnet.model.model import CHGNet
from chgnet.model import StructOptimizer
from pymatgen.core import Structure
from pymatgen.io.cif import *

def run_chgnet(atoms,n_opts=2,rel=StructOptimizer(),relaxer_opts={
	'fmax':[0.1,0.05],
	'steps':[250,750],
	'verbose':[True,True]
	}
	,opt_class=['FIRE','BFGSLineSearch'],
	mode='relax',
	opt_device='cpu',
	use_spglib=True
	):
	
	#First do a quick test to make sure that there are no isolated atoms (> 6 angstroms from their nearest neighbour)
	for i in range(len(atoms)):
		j=list(range(len(atoms)))
		j.remove(i)
		distances=atoms.get_distances(i,j,mic=True)
		if min(distances) > 6.0:
			converged=False
			energy=1.e20
			print(min(distances))
			return atoms,energy,converged
	
	chgnet = CHGNet.load()
	
	write("temp.cif",atoms)
	temp_atoms=Structure.from_file("temp.cif")
	
	if mode == 'relax':
		try:
			for i in range(n_opts):
				relaxer=StructOptimizer(optimizer_class=opt_class[i],use_device=opt_device)
				prediction=relaxer.relax(temp_atoms,fmax=relaxer_opts['fmax'][i],steps=relaxer_opts['steps'][i],verbose=relaxer_opts['verbose'][i])
				temp_atoms=prediction['final_structure']
			
			result=chgnet.predict_structure(temp_atoms)
			energy=result['e']*len(temp_atoms)
			forces=result['f']
			
			temp_atoms.to_file("temp.cif")
			
			if abs(forces.max()) <= relaxer_opts['fmax'][-1]:
				converged=True
			else:
				converged=False
			
			atoms=read("temp.cif")
	   	
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
		
		except:
			converged=False
			energy=1.e20
		
		
	if mode == 'single':
		result=chgnet.predict_structure(temp_atoms)
		energy=result['e']*len(temp_atoms)
		converged=True
	
	return atoms,energy,converged
	
#atoms,energy,converged=run_chgnet(atoms=read("after.cif"))

#print(energy,converged)
