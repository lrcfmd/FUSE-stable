import sys
import os
import glob
import ase
from fuse106.gulp import *
from ase.calculators.gulp import GULP
import platform
import re
from ase.io import *
from threading import Thread
import functools
import time
from func_timeout import func_timeout, FunctionTimedOut
#target = "temp01.res"
#restart= "restart.res"
#head = "head2.txt"

def timeout(timeout):
    def deco(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            res = [Exception('function [%s] timeout [%s seconds] exceeded!' % (func.__name__, timeout))]
            def newFunc():
                try:
                    res[0] = func(*args, **kwargs)
                except Exception as e:
                    res[0] = e
            t = Thread(target=newFunc)
            t.daemon = True
            try:
                t.start()
                t.join(timeout)
            except Exception as je:
                print ('error starting thread')
                raise je
            ret = res[0]
            if isinstance(ret, BaseException):
                raise ret
            return ret
        return wrapper
    return deco


def run_gulp(atoms='',shel=None,kwds='',opts='',lib='',produce_steps='',
	gulp_command='gulp < gulp.gin > gulp.got',gulp_timeout=''):
	iat=len(atoms)
	sleeptime=10
	converged = False
	if os.path.isfile("temp.res"):
		os.remove("temp.res")
	if produce_steps==True:
		try:
			files=glob.glob("atoms*.cif")
			for z in range(len(files)):
				os.remove(files[z])
		except:
			pass
	for i in range(len(kwds)):
		if i ==0:
			if not 'dump temp.res\n' in opts[i]:
				opts[i].append('dump temp.res\n')
			if shel == None:
				calc=GULP(keywords=kwds[i],options=opts[i],library=lib)
			else:
				calc=GULP(keywords=kwds[i],options=opts[i],library=lib,shel=shel)
			atoms.set_calculator(calc)
			atoms.get_calculator().optimized = False
			try:
				if gulp_timeout == '':
					energy=atoms.get_potential_energy()
				
				else:
					#func=timeout(timeout=gulp_timeout)(atoms.get_potential_energy)								
					try:
						energy=func_timeout(gulp_timeout,atoms.get_potential_energy)
						#energy=func()
					
					except:
						time.sleep(sleeptime)
						os.system("taskkill /f /im gulp.exe > kill.txt")
						energy = 1.e20
						converged = False
						#print("timed out on ase gulp call \n\n")
					
				if len(atoms) != iat:
					converged = False
					break
					
				try:
					if glob.glob("gulptmp*") != []:
						if platform.system() == 'Windows':
							files=glob.glob("gulptmp*")
							for z in range(len(files)):
								os.remove(files[z])
						if platform.system() == 'Linux':
							os.system("rm gulptmp*")
				except:
					pass
			except:
				f=open("gulp.gin",'r')
				f2=f.readlines()
				cell=re.findall("\d+\.\d+",f2[6])
				for i in range(len(cell)):
					cell[i]=cell[i][:10]
					
				f.close()
				f3=open("gulp.gin",'w')
				f2[6]=str(str(cell[0])+" "+str(cell[1])+" "+str(cell[2])+" "+str(cell[3])+" "+str(cell[4])+" "+str(cell[5])+"\n")
				for i in range(len(f2)):
					f3.write(f2[i])
				f3.close()
				atoms=read_gulp("gulp.gin")
				atoms.set_calculator(calc)
				try:
					if gulp_timeout == '':
						energy=atoms.get_potential_energy()

					else:
						#func=timeout(timeout=gulp_timeout)(atoms.get_potential_energy)								
						try:
							energy=func_timeout(gulp_timeout,atoms.get_potential_energy)
							#energy=func()
						
						except FunctionTimedOut:
							time.sleep(sleeptime)
							os.system("taskkill /f /im gulp.exe > kill.txt")
							energy = 1.e20
							converged = False
							#print("timed out on ase gulp call \n\n")


					if len(atoms) != iat:
						converged = False
						break
					if glob.glob("gulptmp*") != []:
						if platform.system() == 'Windows':
							os.system("del gulptmp*")
						if platform.system() == 'Linux':
							os.system("rm gulptmp*")
				except:
				#except(ValueError):
					energy=1.0e20
				#	break
			if atoms.get_calculator().optimized == True:
				if atoms.calc.Gnorm <= 0.01:
					#converged = True
					#try:
					#	atoms = read_gulp("temp.res")
					#except:
					#	converged = False
					#	pass
					
					### add in extra check to see if the calculation has REALLY converged
					f4=open("gulp.got",'r').readlines()
					if 'opti' in kwds[i]:
						if '  **** Optimisation achieved ****' in f4:
							converged = True
						else:
							converged = False
			
			if 'sing' in kwds[i]:
				converged = True
			if produce_steps==True:
				label=str("atoms"+str(i+1)+".cif")
				write(label,atoms)
			#os.remove("gulp.got")
			#os.remove("gulp.gin")
		
		if i > 0:
			try:
				f=open("temp.res",'r').readlines()
				new_file=open("gulp.gin",'w')
				for j in f:
					if "title" in j:
						start=f.index(j)-1
					if "totalenergy" in j:
						end=f.index(j)
				
				new_file.write(kwds[i])
				new_file.write("\n")
				for j in opts[i]:
					new_file.write(j)
					new_file.write("\n")
				new_file.write(str("library "+str(lib)))
				new_file.write("\n")
				cp=f[start:end]
				
				for j in cp:
					new_file.write(j)
				new_file.close()
				
				if gulp_timeout == '':
					os.system(gulp_command)
				
				else:
						try:
							energy=func_timeout(gulp_timeout,os.system,args=([gulp_command]))
						
						except FunctionTimedOut:
							time.sleep(sleeptime)
							os.system("taskkill /f /im gulp.exe > kill.txt")
							energy = 1.e20
							converged = False
							#print("timed out on ase gulp call \n\n")

				output=open("gulp.got",'r').readlines()
				if 'opti' in kwds[i]:
					converged = '  **** Optimisation achieved ****\n' in output
				
				if 'sing' in kwds[i]:
					converged = True
					for j in output:
						if "Total lattice energy" in j:
							if "eV" in j:
								e_line=j
					energy = float(e_line[e_line.index('-'):e_line.index('eV')-1])
					atoms = read_gulp("temp.res")
				
				if 'sing' not in kwds[i]:
					try:
						for j in output:
							if "Final energy" in j:
								e_line=j
								break
						
						energy = float(e_line[e_line.index('-'):e_line.index('eV')-1])
						
						atoms=read_gulp("temp.res")
						
					except:
						converged = False
						energy = 1.0e20
			except:
				converged = False
				energy = 1.0e20
				
			
	if converged == True:
		try:
			atoms = read_gulp("temp.res")
		except:
			pass
	else:
		converged = False
		energy = 1.0e20
		
	#view(atoms)
	return atoms, energy, converged

