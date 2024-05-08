from fuse202.all import *
import os

run_fuse(

# FUSE setups ##################################################################
composition={'Ca':3,'Ti':2,'O':7},
max_atoms=50,
imax_atoms=50,
restart=False,
read_exisiting_structures = False, # read in and slice previously generated structures for modules?
path_to_structures = 'reference_structures', # should be a directory containing the structures you want to read in
initial_gen=25,
iterations=50,
check_bonds=True,
dist_cutoff = 1.0,
revert_atoms='',

################################################################################

# Specific bits for the BH search rountine #####################################

################################################################################

search_gen_bh=1,
melt_threshold=150,
rmax=500,
#n_moves={12:1},
T=0.02,

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
gulp_command=os.environ['ASE_GULP_COMMAND']

)
