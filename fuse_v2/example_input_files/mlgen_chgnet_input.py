from fuse202.all import *
import os

run_fuse(

# FUSE setups ##################################################################
composition={'Ca':3,'Ti':2,'O':7},
max_atoms=50,
imax_atoms=50,
restart=False,
read_exisiting_structures = True, # read in and slice previously generated structures for modules?
path_to_structures = 'reference_structures', # should be a directory containing the structures you want to read in
initial_gen=25,
iterations=50,
check_bonds=True,
dist_cutoff = 1.0,
pull_random=False,
pull_spp_rank=True,
use_spglib=True,

################################################################################

# Specific bits for the BH search rountine #####################################

################################################################################

search_gen_bh=1,
melt_threshold=150,
rmax=500,
#n_moves={13:1},

################################################################################

# Specific bits for the Gn-Boss ML model for structure generation ##############

################################################################################


#variables used to run ML structure generation:
#gn boss model for structure generation
generate_gn_boss_structures=True, # if set to true, when FUSE is firt launched, it will run gn-boss to generate the pool of referennce structures for this calculation.
gn_boss_command= os.environ['CONDA']+r"conda run -p"+os.environ['GNBOSS']+"get_cifs_for_FUSE.py", # For my machine, I've setup gn-boss in a seperate python environment, this is the command for that version of python
gn_search='tpe', # 'rand' random search, 'tpe' baysian opt, 'pso' particle swarm
gn_max_step=5000, # number of generation attempts for gn-boss
gn_template_path= os.environ['GNBOSS_TEMP'], #path to template files for using gn-boss 
gn_zn_range=[1,2,3,4], # numbers of formula units to scan with GN-BOSS for generating structures
rank_gn_structures='single', #if None; do not rank structures, this should only be set if pull_random = True above, if 'opti' rank with SPPs AFTER geometry optimising the, if 'sing' rank based on single point calculations with SPPs. 
clear_previous_gn_structures=True, #if set to True, before starting the calcluation, remove any previous structures from reference structures & gn-boss generated_results.
generate_structures_only=False,
ranking='chgnet',

r_kwds=['opti conj conp noelectro','opti conj conp noelectro','opti conp noelectro','opti conp noelectro'], # keywords for the gulp input, 
#multiple strings indicate running GULP in multiple stages
r_gulp_opts=[
['\ninclude ./lib.lib\ndump temp.res\nmaxcyc 250\nstepmax 0.001\ntime 5 minutes\ngmax 0.1\ngtol 0.1\nftol 0.1\nxtol 0.1'],
['\ninclude ./lib.lib\ndump temp.res\nmaxcyc 250\nstepmax 0.5\ntime 5 minutes\ngmax 0.1\ngtol 0.1\nftol 0.1\nxtol 0.1'],
['\ninclude ./lib.lib\ndump temp.res\nmaxcyc 1500\nlbfgs_order 5000\ntime 5 minutes\ngmax 0.1\ngtol 0.1\nftol 0.1\nxtol 0.1'],
['\ninclude ./lib.lib\ndump temp.res\nmaxcyc 1500\nlbfgs_order 5000\ntime 5 minutes\ngmax 0.1\ngtol 0.1\nftol 0.1\nxtol 0.1'],
],	# options for gulp, must include one line per set of inputs in "kwds"
r_lib='dummy.lib', # library file for interatomic potentials				



################################################################################

# Definitions for energy calculator(s) #########################################

################################################################################


### vasp inputs
ctype='mixed', # flag to pass to ase inorder to use gulp as the energy minimiser
calcs=['chgnet','vasp'],
### chgnet inputs
n_opts=3,
rel=StructOptimizer(),
relaxer_opts={
'fmax':[0.5,0.1,0.05],
'steps':[500,500,1000],
'verbose':[True,True,True]
},
opt_class=['SciPyFminCG','FIRE','BFGSLineSearch'],
opt_device='cpu',


vasp_opts={
'1':Vasp(prec='Fast',ediff=1.E-4,nelm=50,encut=600,ibrion=2,isif=3,nsw=100,ediffg=0.05,nwrite=1,ncore=10,algo='Fast',gamma=True,lreal='Auto',setups={'Li':'_sv','Sn':'_d'},kspacing=0.5),
'2':Vasp(prec='Fast',ediff=1.E-4,nelm=50,encut=600,ibrion=2,isif=3,nsw=100,ediffg=0.05,nwrite=1,ncore=10,algo='Fast',gamma=True,lreal='Auto',setups={'Li':'_sv','Sn':'_d'},kspacing=0.5),
'3':Vasp(prec='Normal',ediff=1.E-5,nelm=50,encut=650,ibrion=1,isif=3,nsw=20,ediffg=-0.03,nwrite=1,ncore=10,gamma=True,lreal='.False.',setups={'Li':'_sv','Sn':'_d'},kspacing=0.3,ismear=0),
'4':Vasp(prec='Normal',ediff=1.E-5,nelm=50,encut=650,ibrion=1,isif=3,nsw=20,ediffg=-0.03,nwrite=1,ncore=10,gamma=True,lreal='.False.',setups={'Li':'_sv','Sn':'_d'},kspacing=0.3,ismear=0),
'5':Vasp(prec='Normal',ediff=1.E-5,nelm=50,encut=650,ibrion=1,isif=3,nsw=150,ediffg=-0.03,nwrite=1,ncore=10,gamma=True,lreal='.False.',setups={'Li':'_sv','Sn':'_d'},kspacing=0.3,ismear=0),
},
kcut=[20],

lib='dummy.lib', # library file for interatomic potentials
shel=[''],      # species using shells in gulp
gulp_timeout='', # timout command for running gulp, leave as '' if not using, only use in Windows!
assemble_spp_=True, # If set to true collate SPP potential for the system
spp_path=os.environ['SPP_PATH'],
gulp_command=os.environ['ASE_GULP_COMMAND'],


)
