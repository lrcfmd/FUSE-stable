import os
import ase
from ase import *
from ase.io import *
import spglib
import sys
from ase.calculators.castep import Castep

def run_castep(atoms='',castep_opts='',dist_cutoff=1.0,use_spglib=True,ftol=0.5,pspot_setting=''):

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

	if short_contact == False:
		new_atoms=atoms.copy()
		for i in range(len(list(castep_opts.keys()))):
			#print(new_atoms)
			
			calc=castep_opts[str(list(castep_opts.keys())[i])]
			#calc._find_pspots=True
			#calc.param.write_cell_structure=True
			#calc._rename_existing_dir=False
			#calc._build_missing_pspots=True
			#calc._copy_pspots=True
			##path=str(os.getcwd()+'/CASTEP/')
			#calc._castep_pp_path=os.environ['CASTEP_PP_PATH']
			
			if pspot_setting != '':
				calc.find_pspots(pspot=pspot_setting)
				#for x in range(len(pspot_setting['elements'])):
				#	#lab=pspot_setting['elements'][x]+"_"+pspot_setting['type']+"_"+pspot_setting['xc']+"_"+pspot_setting['suff']+"."+pspot_setting['extension']
				#	calc.cell.species_pot=(pspot_setting['elements'][x],pspot_setting['type'][x])
					
			#print(ps)	
			#print(calc)
			
			if os.path.isfile('CASTEP/castep.castep'):
				os.remove('CASTEP/castep.castep')

			if use_spglib == True:
				try:    
					lattice,positions,numbers=spglib.standardize_cell(new_atoms,symprec=1.e-5)
					temp2=Atoms(numbers=numbers,pbc=True)
					temp2.cell=lattice
					temp2.set_scaled_positions(positions)
					new_atoms=temp2.copy()
					#print("I'm using SPGLIB!")
				except:
					#print("I failed at using SPGLIB!!")
					pass

			new_atoms.calc = calc
						
			try:
				if new_atoms.calc.dryrun_ok():
					dry_ok=True
					
			except:
				print('failed dry run, castep doesnt like something')
				converged = False
				energy = 1.e20

			#sys.exit()

			if dry_ok==True:
				#print("\n",dir(castep_opts[str(list(castep_opts.keys())[i])].param.task._value))
				#print("hello ", castep_opts[str(list(castep_opts.keys())[i])].param.task._value=='SinglePoint')
				#break
				
				try:
					
					energy = new_atoms.get_potential_energy()
					#print("\n\nI've run castep ok, energy = ",str(energy),"\n\n")
					new_atoms=read('CASTEP/castep-out.cell')
					
					with open('CASTEP/castep.castep', 'r') as out:
						lines = out.readlines()
						for i in lines:
							if 'Geometry optimization completed successfully' in i:
								converged=True
								break
							else:
								converged=False
								

				except:
					if str(calc.param.task._value) != 'SinglePoint':
						converged = False
						with open('CASTEP/castep.castep','r') as out:
							lines=out.readlines()
							for i in lines:
								if 'finished' in i:
									energy=float(str(i)[48:64])
									new_atoms = read('CASTEP/castep.geom')
									break
								else:
									converged=False
									energy=1.e20
					else:
						print('failed single point calculation, something is very wrong')
						converged = False
						energy = 1.e20
						
			if str(calc.param.task._value) == "SinglePoint":
				converged=True
						
		

		temp_atoms=new_atoms.repeat([2,2,2])
		temp1=temp_atoms.get_all_distances()
		temp2=[]
		for i in range(len(temp1)):
			for j in range(len(temp1[i])):
				if temp1[i][j] != 0:
					temp2.append(temp1[i][j])
		if min(temp2) <= dist_cutoff:
			converged=False

	        # also put a check in for forces if the calculation isn't converged
		if str(calc.param.task) != 'SinglePoint':
			if converged == False:
				try:
					with open('CASTEP/castep.castep','r') as output:
						lines=output.readlines()
						for i in lines:
							if '|F|max' in i:
								fo=float(i[16:30])
					print("max forces: ",fo)
					if fo > ftol:
						converged = False
						energy = 1.e20
				except:
					converged = False
					energy = 1.e20
		#print(new_atoms, energy, converged)
	return new_atoms, energy, converged

#atoms=read('castep.cell')
#print(atoms)
#castep_opts={
#'1':Castep(task='GeometryOptimisation',kpoint_mp_spacing=0.5,xc_functional='PBE',basis_precision='coarse',max_scf_cycles=14,mix_charge_amp=0.3,
#geom_energy_tol=1,geom_force_tol=2,geom_stress_tol=20,geom_disp_tol=0.003)
#,'2':Castep(task='SinglePoint',kpoint_mp_spacing=0.5,write_cell_structure=True,elec_energy_tol=0.0005,xc_functional='PBE',cut_off_energy=500)
#}
#res = run_castep(atoms=atoms,castep_opts=castep_opts)
#print(res[0],res[1],res[2])

