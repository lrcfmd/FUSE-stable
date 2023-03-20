from fuse106 import *
import sys
import time

#### composition information
#composition = {'Y':[2],'Ti':[2],'O':[7]} #composition in non-charged format
composition = {'Y':[2,+3],'Ti':[2,+4],'O':[7,-2]} #composition in ionic format
#composition = {'Sr':[1],'Ti':[1],'O':[3]}
#composition = {'Sr':[1,+2],'Ti':[1,+4],'O':[3,-2]}

#### set search routine 1 = BH, 2 = GA
search = 1

### Genetic algo inputs (running on defaults)

### Basin hopping inputs (running on defaults)
#n_moves={1:16,2:21,3:1,4:21,5:10,6:1,7:1,9:16,10:16,11:10}
n_moves={1:1,2:1,3:1,4:1,5:1,6:1,7:1,9:1,10:1,11:1}
#n_moves={11:1}
r_moves={6:2,7:1,10:4}
#n_moves={7:1}
### general inputs
rmax=7000 # the number of steps since the global minimum was located before ending
#the calculation
iterations=250 # number of structures to compute in this run
restart=False # restart a previous run?
initial_gen=100 # size of the initial population
search_gen=1 # size of each generation in to be used by the search routine
max_atoms=50 # global limit on the maximum number of atoms
imax_atoms=150 # maximum number of atoms to be used in the initial population
check_bonds=True
check_distances=True
melt_threshold=0.15*rmax

### gulp inputs
ctype='gulp' # flag to pass to ase inorder to use gulp as the energy minimiser
kwds=['opti conp conj c6','opti conp c6'] # keywords for the gulp input, 
#kwds=['sing conj conp c6'] # keywords for the gulp input, 
#multiple strings indicate running GULP in multiple stages
gulp_opts=[
['\nlibrary lib2.lib\ndump temp.res\ntime 15 minuets\nmaxcyc 75\nstepmax 0.075\ntime 5 minutes'],
['\nlibrary lib2.lib\ndump temp.res\ntime 15 minutes\nmaxcyc 2000\ntime 10 minutes\ngmax 0.1\ngtol 0.1\nxtol 0.1\nftol 0.1'],
]	# options for gulp, must include one line per set of inputs in "kwds"
lib='lib2.lib' # library file for interatomic potentials
shel=['']	# species using shells in gulp
gulp_timeout=500 # timeout for gulp calls in seconds, leave as '' to ignore, use on Windows machines only.

# output options (running on defaults)
write_graph=True
write_all_structures=True
# run fuse

run_fuse(composition=composition,search=search,rmax=rmax,
	iterations=iterations,restart=restart,initial_gen=initial_gen,
	search_gen=search_gen,max_atoms=max_atoms,imax_atoms=imax_atoms,ctype=ctype,
	kwds=kwds,gulp_opts=gulp_opts,lib=lib,shel=shel,gulp_timeout=gulp_timeout,write_graph=write_graph,check_bonds=check_bonds,write_all_structures=write_all_structures,
	check_distances=check_distances,n_moves=n_moves,
	melt_threshold=melt_threshold,r_moves=r_moves
	)

sys.exit()
