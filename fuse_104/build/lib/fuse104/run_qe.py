from ase import *
from ase.io import read, write
from ase.calculators.espresso import Espresso
import math
import os
import sys
import glob

##########################################################################################
import shlex
import re
import numpy
from numpy import arccos, pi, dot
from numpy.linalg import norm

def cellpar(atoms):
	cell = atoms.cell
	a = norm(cell[0])
	b = norm(cell[1])
	c = norm(cell[2])
	alpha = arccos(dot(cell[1], cell[2])/(b*c))*180./pi
	beta  = arccos(dot(cell[0], cell[2])/(a*c))*180./pi
	gamma = arccos(dot(cell[0], cell[1])/(a*b))*180./pi

	cell = []
	cell=[a,b,c,alpha,beta,gamma]
	return cell
##########################################################################################


def check_convergence(output_file=''):
    converged = False
    lines = open(output_file, 'r').readlines()
    if '   JOB DONE.\n' not in lines:
        converged = False
        return converged
    if ('nstep' in j for j in lines):
        if 'A final scf calculation at the relaxed structure.' in lines:
            converged = True
    elif ('nstep' not in j for j in lines):
        if ('convergence has been achieved' in j for j in lines):
            converged = True
    return converged


def run_qe(atoms='', qe_opts='', kcut='', produce_steps=''):

    new_atoms=atoms.copy()
    for i in range(len(list(qe_opts.keys()))):
        if len(kcut) >= 2:
            temp_kcut=kcut[i]
        else:
            temp_kcut=kcut
        cell=new_atoms.get_cell_lengths_and_angles()
        kp=[]
        kp=(int(math.ceil(temp_kcut/cell[0])),
            int(math.ceil(temp_kcut/cell[1])),
            int(math.ceil(temp_kcut/cell[2])))
        calc=qe_opts[str(list(qe_opts.keys())[i])]
        calc.set(kpts=kp)
        calc=qe_opts[str(list(qe_opts.keys())[i])]
        new_atoms.set_calculator(calc)
        try:
            energy=new_atoms.get_potential_energy()
        except:
            converged = False
            energy = 1.e20
        if produce_steps==True:
            label=str('atoms'+str(i+1)+'.cif')
            write(label,new_atoms)
    
        # remove any output / tempory files before starting"
        old_files=glob.glob("pwscf*")
        for k in range(len(old_files)):
            os.system("rm -r "+str(old_files[k]))
        
    converged=check_convergence('espresso.pwo')
#    print("converged:" +str(converged)+str("\n\n"))
	
    temp_atoms=new_atoms.repeat([2,2,2])
    temp1=temp_atoms.get_all_distances()
    temp2=[]
    for i in range(len(temp1)):
        for j in range(len(temp1[i])):
            if temp1[i][j] != 0:
                temp2.append(temp1[i][j])
    if min(temp2) <= 1.2:
        converged=False

    energy = energy/13.6056980659 #convert from Ry to eV
    return new_atoms,energy,converged	


#qe_opts={
#'1':Espresso(pseudopotentials={'Sr':'Sr.pbe-spn-kjpaw_psl.1.0.0.UPF',
#                               'Nb':'Nb.pbe-spn-kjpaw_psl.1.0.0.UPF',
#                               'O':'O.pbe-n-kjpaw_psl.1.0.0.UPF'},
#             input_data={'control':{'calculation':'vc-relax',
#                                    'pseudo_dir':'./pp/',
#                                    'outdir':'./',
#                                    'etot_conv_thr':0.0001,
#                                    'forc_conv_thr':0.001,
#                                    'disk_io':'low',
#                                    'nstep':50},
#                         'system':{'ecutwfc':29.39945,
#                                   'occupations':'smearing',
#                                   'smearing':'gaussian',
#                                   'degauss':0.0073498},
#                         'electrons':{'electron_maxstep':50,
#                                      'conv_thr':0.000001}},
#             koffset=(0,0,0))
#}
#
#atoms,energy,converged=run_qe(atoms=read('SrNbO3_converged.cif'),qe_opts=qe_opts,kcut=[20])
#
#print(str('energy: ')+str(energy))
#print(str('converged: ')+str(converged))
