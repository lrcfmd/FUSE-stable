from fuse104.all import *
import sys
import time

#### composition information
#composition = {'Y':[2],'Ti':[2],'O':[7]} #composition in non-charged format
composition = {'Sr':[1,+2],'Ti':[1,+4],'O':[3,-2]} #composition in ionic format

search = 1 # basin hopping search


### general inputs
rmax=10000 # the number of steps since the global minimum was located before ending
#the calculation
iterations=1 # number of structures to compute in this run
restart=False # restart a previous run?
initial_gen=1 # size of the initial population
search_gen=1 # size of each generation in to be used by the search routine
max_atoms=5 # global limit on the maximum number of atoms
imax_atoms=5 # maximum number of atoms to be used in the initial population
check_bonds=False

### gulp calculator options
ctype='mixed' # flag to pass to ase inorder to use gulp as the energy minimiser
calcs=['gulp','vasp'] # which calculators to use and in which order

# gulp options
kwds=['opti conj conp c6','opti lbfgs conp c6'] # keywords for the gulp input, 
#multiple strings indicate running GULP in multiple stages
gulp_opts=[
['\nlibrary lib2.lib\ndump temp.res\ntime 15 minuets\nmaxcyc 50\nstepmax 0.1\ntime 5 minutes'],
['\nlibrary lib2.lib\ndump temp.res\ntime 15 minutes\nmaxcyc 1500\nlbfgs_order 5000\ntime 5 minutes'],
]	# options for gulp, must include one line per set of inputs in "kwds"
lib='lib2.lib' # library file for interatomic potentials
shel=['']	# species using shells in gulp

#vasp options
vasp_opts={
'1':Vasp(prec='Fast',ediff=1.E-4,nelm=50,encut=400,ibrion=2,isif=3,nsw=50,ediffg=0.05,nwrite=1,algo='Fast',gamma=True,lreal='Auto',setups={'Sr':'_sv'}),
'2':Vasp(prec='Fast',ediff=1.E-4,nelm=50,encut=400,ibrion=2,isif=3,nsw=50,ediffg=0.05,nwrite=1,algo='Fast',gamma=True,lreal='Auto',setups={'Sr':'_sv'}),
'3':Vasp(prec='Normal',ediff=1.E-4,nelm=50,encut=450,ibrion=1,isif=3,nsw=250,ediffg=-0.05,nwrite=1,gamma=True,lreal='Auto',setups={'Sr':'_sv'}),
}
kcut=[20]
# output options (running on defaults)

# run fuse

run_fuse(composition=composition,search=search,rmax=rmax,
	iterations=iterations,restart=restart,initial_gen=initial_gen,
	search_gen=search_gen,max_atoms=max_atoms,imax_atoms=imax_atoms,ctype=ctype,
	kwds=kwds,gulp_opts=gulp_opts,lib=lib,shel=shel,calcs=calcs,check_bonds=check_bonds,
        vasp_opts=vasp_opts,kcut=kcut)

sys.exit()
