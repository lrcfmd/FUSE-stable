from fuse202.all import *
import os

run_fuse(

# FUSE setups ##################################################################
composition={'Ca':3,'Ti':2,'O':7}, # the emperical formula unit for this calculation.
max_atoms=50,  # Maximum number of atoms to use in the overall calculation
imax_atoms=50, # Maximum number of atoms to use in only the initial population, for large systems, it can be helpful to set this smaller than the "max_atoms".
restart=False,  # restart a previous calculation? If set to True, FUSE will attempt to restart from a previous calculation.
read_exisiting_structures = False, # read in and slice previously generated structures for modules?
path_to_structures = 'reference_structures', # should be a directory containing the structures you want to read in
initial_gen=25, # the number of structures to include in the initial population
iterations=15,  # the number of structures to perform energy calculations on in this run of FUSE.

################################################################################

# Specific bits for the BH search rountine #####################################

################################################################################

search_gen_bh=1, # Number of new structures to generate at each step in the basin hopping search
melt_threshold=150, # After this number of steps since the bashin hopping routine makes a downhill step, FUSE will increase the temperateure parameter to escape local minima
rmax=500, # The basin hopping routine will be considered converged if this many structures are generated since the current lowest energy structure was located.
T=0.02, # the temperature parameter for the basin hopping routine, the higher the value, the larger uphill step that the basin hopping routine will accept.

################################################################################

# Definitions for energy calculator(s) #########################################

################################################################################

### gulp inputs
ctype='mixed', # flag to pass to ase inorder to use gulp as the energy minimiser

calcs=['gulp','vasp'],

# Options for VASP

# For more information to configure Vasp, please see the ase documentation available here:
# https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#vasp

vasp_opts={
'1':Vasp(prec='Fast',ediff=1.E-4,nelm=50,encut=520,ibrion=2,isif=3,nsw=100,ediffg=0.05,nwrite=1,ncore=18,algo='Fast',gamma=True,lreal='Auto',setups={'Ca':'_sv'},kspacing=0.5),
'2':Vasp(prec='Fast',ediff=1.E-4,nelm=50,encut=520,ibrion=2,isif=3,nsw=50,ediffg=0.05,nwrite=1,ncore=18,algo='Fast',gamma=True,lreal='Auto',setups={'Ca':'_sv'},kspacing=0.5),
'3':Vasp(prec='Normal',ediff=1.E-5,nelm=50,encut=520,ibrion=1,isif=3,nsw=10,ediffg=-0.03,nwrite=1,ncore=18,gamma=True,lreal='.False.',setups={'Ca':'_sv'},kspacing=0.3,ismear=0),
'4':Vasp(prec='Normal',ediff=1.E-5,nelm=50,encut=520,ibrion=1,isif=3,nsw=10,ediffg=-0.03,nwrite=1,ncore=18,gamma=True,lreal='.False.',setups={'Ca':'_sv'},kspacing=0.3,ismear=0),
'5':Vasp(prec='Normal',ediff=1.E-5,nelm=50,encut=520,ibrion=1,isif=3,nsw=50,ediffg=-0.03,nwrite=1,ncore=18,gamma=True,lreal='.False.',setups={'Ca':'_sv'},kspacing=0.3,ismear=0),
},
kcut=[20], #This is now a redundent parameter, please leave it here as a dummy value. this will be removed in a later version.

#Options for GULP
# For more details on how to use GULP with ASE please see the documentation here:
# https://wiki.fysik.dtu.dk/ase/ase/calculators/gulp.html

kwds=['opti conp conj noelectro','opti conp noelectro conj','opti conp conj noelectro','opti conp conj noelectro','sing conp noelectro'], # keywords for the gulp input,

gulp_opts=[
['\ndump temp.res\nmaxcyc 250\ntime 5 minutes\ncutp 11.0 2.0\ninclude ./lib.lib\nstepmax 0.001'],
['\ndump temp.res\nmaxcyc 250\ntime 5 minutes\ncutp 11.0 2.0\ninclude ./lib.lib\ngmax 0.1\ngtol 0.1\nftol 0.1\nxtol 0.1\stepmax 0.5'],
['\ndump temp.res\nmaxcyc 1500\ntime 5 minutes\ncutp 11.0 2.0\ninclude ./lib.lib\ngmax 0.1\ngtol 0.1\nftol 0.1\nxtol 0.1'],
['\ndump temp.res\nmaxcyc 1500\ntime 5 minutes\ncutp 11.0 2.0\ninclude ./lib.lib\ngmax 0.1\ngtol 0.1\nftol 0.1\nxtol 0.1\nswitch_minimiser bfgs gnorm 1'],
['\ndump temp.res\nmaxcyc 1500\ntime 5 minutes\ncutp 11.0 2.0\ninclude ./lib.lib\ngmax 0.1\ngtol 0.1\nftol 0.1\nxtol 0.1\nswitch_minimiser bfgs gnorm 1'],

],   # options for gulp, must include one line per set of inputs in "kwds"
lib='dummy.lib', # library file for interatomic potentials
shel=[''],   # species using shells in gulp
# GULP options for FUSE:

gulp_timeout='', # timeout for gulp calls in seconds, leave as '' to ignore, use on Windows machines only.
assemble_spp_=True, # create a fresh SPP potential for this calculation, for most applications, this should be set to True.
spp_path=os.environ['SPP_PATH'], # if you have set teh enviroment variables required in the README file for FUSE v.2, you should not need to change this.
gulp_command=os.environ['ASE_GULP_COMMAND'], # if you have set teh enviroment variables required in the README file for FUSE v.2, you should not need to change this.
)
