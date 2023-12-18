# function to take the updated structure object and return a new atoms object
# this should be used after a structure has been altered

from fuse202.all import *

# example to take an input cif and break it down into modules to test the re-assembly
##name="garnet"
##
##bondtable=numpy.load("bondtable.npz",allow_pickle=True)
##bondtable=bondtable['bond_table'].item()
##
##input_files=[ 'Ca3Al2Si3O12.cif',
##             #'Sr4Ti3O10_34630.cif'
##             #'Y2O3_23811.cif'
##            ]
##
##structure=extract_module(input_files,bondtable)
##
##print(structure.keys())


#print(cell_par)

#this is assuming that you have both the "structure" object for the structure that you want to assemble
def assemble_structure2(structure):
	modules=structure['modules']
	sub_mod_cell=structure['sub module cell']
	cell_shape=structure['shape in submods']
	ap=structure['ap']
	# 1. set the target cell shape in sub-modules
	cell_par=[round(ap*cell_shape[0],6),round(ap*cell_shape[1],6),round((ap/2)*cell_shape[2],6),cell_shape[3],cell_shape[4],cell_shape[5]]
	#view(modules)
	
	
	# 2. start pulling out modules
	
	position=[0,0,0]
	
	#lets start by step by step building the unit cell
	
	atoms=Atoms()
	atoms.pbc=['True','True','True']
	try:
		atoms.cell=cell_par
	except:
		atoms.cell=[cell_par[0],cell_par[1],cell_par[2],90,90,90]
	#print(atoms.cell.cellpar())
	#view(atoms)
	
	nsub=0
	
	for i in modules:
	    #print("position: ", position)
	    translation=[ position[0]*ap, position[1]*ap, position[2]*(ap/2) ]
	    #print("translation: ", translation)
	    #print(i.get_positions())
	    i.translate(translation)
	    #print(i.get_positions())
	    atoms = atoms+i
	    
	    nsub+=1
	    #write(str("layer_"+str(nsub)+".cif"),atoms)
	    if nsub < len(modules):
	        
	        if position[0] < cell_shape[0]-1:
	            position[0] += 1
	            continue
	        
	        if position[0] == cell_shape[0]-1:
	            
	            if position[1] == cell_shape[1]-1:
	                position[2] += 1
	                position[0] = 0
	                position[1] = 0
	                continue
	                
	            if position[1] < cell_shape[1]-1:
	                position[1] += 1
	                position[0] = 0
	return atoms              


# In[ ]:

# test code to take the assembled structure and optimise it with gulp

#write(str("assembled_"+str(name)+".cif"),atoms) 
#
#from ase.calculators.gulp import GULP
#
#label=str('output cif relaxed_'+str(name)+"_REF")
#calc=GULP(keywords='opti conp',options=[label],library='pedone.lib')
#
#ref_atoms=read(input_files[0])
#ref_atoms.set_calculator(calc)
#ref_energy=ref_atoms.get_potential_energy()
#print(ref_energy/len(ref_atoms))
#
#label=str('output cif relaxed_'+str(name))      
#calc=GULP(keywords='opti conp',options=[label],library='pedone.lib')
#atoms.set_calculator(calc)
#energy=atoms.get_potential_energy()
#print(energy/len(atoms))
#
#
# In[ ]:




