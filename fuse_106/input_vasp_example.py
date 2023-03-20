from fuse106.all import *
import sys
import time

#### composition information
composition = {'Sr':[1,+2],'Ti':[1,+4],'O':[3,-2]} #composition in ionic format

#### set search routine 1 = BH, 2 = GA, 3 = ML
search = 1

### Genetic algo inputs

### Basin hopping inputs

### general inputs
rmax=150
iterations=2
restart=False
initial_gen=1
search_gen=1
max_atoms=10
imax_atoms=10

### vasp inputs
ctype='vasp'
vasp_opts={
'1':Vasp(prec='Fast',ediff=1.E-4,nelm=50,encut=400,ibrion=2,isif=3,nsw=50,ediffg=0.05,nwrite=1,ncore=12,algo='Fast',gamma=True,lreal='Auto',setups={'Sr':'_sv'}),
'2':Vasp(prec='Fast',ediff=1.E-4,nelm=50,encut=400,ibrion=2,isif=3,nsw=50,ediffg=0.05,nwrite=1,ncore=12,algo='Fast',gamma=True,lreal='Auto',setups={'Sr':'_sv'}),
'3':Vasp(prec='Normal',ediff=1.E-4,nelm=50,encut=450,ibrion=1,isif=3,nsw=250,ediffg=-0.05,nwrite=1,ncore=12,gamma=True,lreal='Auto',setups={'Sr':'_sv'}),
}
kcut=[20]

# output options

# run fuse
run_fuse(composition=composition,search=search,initial_gen=initial_gen,
	  max_atoms=max_atoms,ctype=ctype,vasp_opts=vasp_opts,kcut=kcut,
	  iterations=iterations,restart=restart,search_gen=search_gen,rmax=rmax,
     imax_atoms=imax_atoms)

sys.exit()
