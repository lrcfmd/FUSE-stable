from fuse205.all import *

run_fuse(

# FUSE setups ##################################################################
composition={'Sr':1,'Ti':1,'O':3}, # the emperical foruma unit for this calculation
max_atoms=10, # Maximum number of atoms to use in the overall calculation
imax_atoms=10, # Maximum numbe rof atoms to use in only the inital population, for large systems,it can be helpful to set this smaller than the "max_atoms".
restart=False, # restart a previous calculation? If set to True, FUSE will attempt to restart from a previous calculation.
read_exisiting_structures = True, # read in and slice previously generated structures for modules?
path_to_structures = 'reference_structures', # should be a directory containing the structures you want to read in
initial_gen=1, # the number of structures to include in the initial population
iterations=10, # the number of structures to perform energy calculations on in this run of FUSE.
pull_random=False, # if set to True, when using pre-built structure pull them in a random order
pull_spp_rank=True, # if True, pull pre-built structures as ranked by spp potentials. BEWARE! you need to first run a script to rank all pre-generated structures & produce one or more csv files containing file names and corresponding energies
use_spglib=True, # if True, everytime an atoms object is generated, tidy up the structure with spglib before the modules are extracted 
search=2, # # search routine for FUSE to use, 1 = basin hopping, 2 basin hopping with reinforcement learning

################################################################################

# Specific bits for the BH search rountine #####################################

################################################################################

search_gen_bh=1, # Number of new structures to generate at each step in the basin hopping search
melt_threshold=150, # After this number of steps since the bashin hopping routine makes a downhill step, FUSE will increase the temperateure parameter to escape local minima
rmax=1500, # The basin hopping routine will be considered converged if this many structures are generated since the current lowest energy structure was located. 

################################################################################

# Specific bits for structure generation ##############

################################################################################

generator=['gnboss'], # List of algorithms to use for structure generation
rank_structures='single', # if None; do not rank structures, this should only be set if pull_random = True above, if 'opti' rank with SPPs AFTER geometry optimising the, if 'sing' rank based on single point calculations with SPPs.
generate_structures_only=False, # Run the structure generation & ranking, then exit if set to True. To then re-start with the basin hopping routine, set restart above to False along with generate_gn_boss_structures and/or generate_airss_structures and clear_previous_structures. 
ranking='chgnet', # if None; do not rank structures, this should only be set if pull_random = True above, if 'opti' rank AFTER geometry optimising the, if 'sing' rank based on single point calculations.
clear_previous_structures=True, # if set to True, before starting the calcluation, remove any previous structures from reference structures, gn-boss generated results and airss generate results.

# Specific bits for GN-Boss model for structure generation
generate_gn_boss_structures=True, # if set to true, when FUSE is firt launched, it will run gn-boss to generate the pool of referennce structures for this calculation.
gn_boss_command= os.environ['GNBOSS']+'/bin/python get_cifs_for_FUSE.py', # For my machine, I've setup gn-boss in a seperate python environment, this is the command for that version of python
gn_search='rand', # 'rand' random search, 'tpe' baysian opt, 'pso' particle swarm
gn_max_step=5000, # number of generation attempts for gn-boss
gn_template_path= os.environ['GNBOSS_TEMP'], #path to template files for using gn-boss 
gn_zn_range=[1], # numbers of formula units to scan with GN-BOSS for generating structures

################################################################################

# Definitions for energy calculator(s) #########################################

################################################################################

ctype='chgnet', 

### chgnet inputs
n_opts=3,
rel=None,
relaxer_opts={
'fmax':[0.5,0.1,0.05],
'steps':[500,500,1000],
'verbose':[True,True,True]
},
opt_class=['FIRE','FIRE','BFGSLineSearch'],
opt_device='cpu',

assemble_spp_=False, # If set to true collate SPP potential for the system
spp_path=os.environ['SPP_PATH'],

######### Specific bits for reinforce algorithm #############

#params_db={'host': 'localhost','database': 'testing_db','user': 'root','password': ''}, # uncomment this to use mysql database for storing reinforce data 
params_db={'database': 'TEST.db'}, # uncomment this if you want to store reinforce data using sqlite3. Setting database as a .db filename creates a file in the working directory, setting database as ':MEMORY:' doesnt and stores the data elsewhere. 
reinforce_verbosity=False, # set as True to print all output text to do with reinforce, there's a lot!
target_energy='', # set a target energy that will end the run if met, for testing purposes

### mysql specific options
dir_num = '1', # numerical identifier for mysql tables
reinforce_table = 'test_table_', # name for the table, set inputs as the same across runs to cross_learn if using mysql
reinforce_theta_table = 'test_table_theta_', # as above but for theta

alpha=0.0005,
reg_params={'e_threshold': 0.8, 'beta': 1, 'free_term': 1,
	'h_type': 'linear', 'zero_reward_penalty': -3,
	'non_converge_penalty': -3,
	'non_unique_penalty': -3,
	'scale_reward': True,
        'epsilon': 0.1,
	'non_unique_reward': False,
        'step_reward_limit': 5000,
        'reward_limit_last': True,
        'smart_penalty': True,
	'storage_type': 'sqlite'} # 'mysql' if using mysql to store data, 'sqlite' if using sqlite

)

