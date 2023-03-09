from fuse_stable import *
import sys
import time

#### composition information
composition = {'Sr':[1,+2],'Ti':[1,+4],'O':[3,-2]} #Setting the composition for the calculation

### general inputs
rmax=50 # the number of steps since the global minimum was located before ending
#the calculation
iterations=10 # number of structures to compute in this run
restart=False # restart a previous run?
initial_gen=5 # size of the initial population
search_gen=1 # size of each generation in to be used by the search routine
max_atoms=20 # global limit on the maximum number of atoms
imax_atoms=20 # maximum number of atoms to be used in the initial population

### gulp inputs ### you shouldn't need to edit below this line!
ctype='gulp' # flag to pass to ase inorder to use gulp as the energy minimiser
kwds=['opti conj conp c6','opti conp c6'] # keywords for the gulp input, 
#multiple strings indicate running GULP in multiple stages
gulp_opts=[
['\nlibrary lib2.lib\ndump temp.res\nmaxcyc 50\nstepmax 0.1\ntime 5 minutes'],
['\nlibrary lib2.lib\ndump temp.res\nmaxcyc 1500\nlbfgs_order 5000\ntime 5 minutes'],
]	# options for gulp, must include one line per set of inputs in "kwds"
lib='lib2.lib' # library file for interatomic potentials
shel=['']	# species using shells in gulp

# output options (running on defaults)

# run fuse

run_fuse(composition=composition,rmax=rmax,
	iterations=iterations,restart=restart,initial_gen=initial_gen,
	search_gen=search_gen,max_atoms=max_atoms,imax_atoms=imax_atoms,ctype=ctype,
	kwds=kwds,gulp_opts=gulp_opts,lib=lib,shel=shel)

sys.exit()
