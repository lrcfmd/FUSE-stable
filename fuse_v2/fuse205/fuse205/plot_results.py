import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas
import sys
import numpy as np
import glob
import os
from ase.io import *

def plot_graph(search=1,initial_population=''):
#graph_output={'move':[],'type':[],'step':[],'energy':[],'temperature':[],'current energy':[],'global minimum energy':[],'structure file name':[],'accepted?':[] }		
	#create sactter plot
	data=pandas.read_csv("graph_output.csv")
	#print (data)
	n=0
	for i in range(len(data['step'])):
		if data['current energy'][i] == 0:
			n+=1
	
	fig1, ax1=plt.subplots()
	colour='tab:red'
	ax1.set_xlabel('Structures')
	ax1.set_ylabel('Energy (eV/atom)',color='tab:red')
	ax1.scatter(data['step'],data['energy'],color=colour,s=0.75)
	ax1.plot(data['step'][n:],data['current energy'][n:],color='tab:green',linewidth=1)
	ax1.plot(data['step'][n:],data['global minimum energy'][n:],color='m',linewidth=0.75)
	lim=min(data['energy']) +max(data['energy'])-min(data['energy'])
	
	if abs( min(data['energy']) +max(data['energy'])-min(data['energy']) )  >= 0.5:
		lim=min(data['energy'])+0.5
	
	#print("data limit: ",lim)
	
	ax1.set_ylim([ min(data['energy'])-0.025, lim ])
	ax1.tick_params(axis='y',labelcolor=colour)
	#ax1.legend(loc=2,framealpha=1)
	#colour='tab:green'
	#ax3=ax1.twinx()
	#ax3.set_ylabel('Energy (meV/atom)')
	#ax3.tick_params(axis='y',labelcolor='tab:red')
	if search == 1:
		ax2=ax1.twinx()
		#colour='tab:black'
		ax2.set_ylabel("Theta (eV)",color='k')
		ax2.plot(data['step'],data['temperature'],color='k',linewidth=0.75)
		ax2.tick_params(axis='y',labelcolor='k')	
		
	fig1.tight_layout()
	#plt.show()
	plt.savefig('run_output.png',dpi=600)
	del fig1,ax1
	try:
		del ax2
	except:
		pass
	
	#create plot of used moves distribution
	used_moves={1:0,2:0,3:0,4:0,5:0,6:0,7:0,8:0,9:0,10:0,11:0,12:0,13:0,14:0}
	
	for i in data['move']:
		if not i == "initial":
			used_moves[int(i)]+=1
	
	used_moves2={}
	for i in list(used_moves.keys()):
		if used_moves[i] > 0:
			used_moves2[i]=used_moves[i]
	
	labels=list(used_moves2.keys())
	values=list(used_moves2.values())
	plt.subplots()
	plt.pie(values,labels=labels)
	plt.savefig("moves_used.png",dpi=600)
	#plt.show()
	
	
	
	#create plot of all structures, energy vs. density
	data3={'density':[],'energy':[]}
	#files=data['structure file name']
	#files=glob.glob("../structures/*.cif")
	energies=data['energy']
	#print("\n\nhello!")
	#print("\n\n",os.getcwd())
	#temp_dat=pickle.load(open("../restart/initial_population.p",'rb'))
	temp_dat=initial_population
	#print("hello")
	keys=list(temp_dat.keys())
	for i in keys:
		eng=temp_dat[i]['energy']
		if eng > 20:
			continue
		#f=files[i]
		#path="../structures/"+f
		at=temp_dat[i]['atoms'].copy()
		
		amu_to_g=1.6605E-24
		ang3_to_cm3=1e-24
		
		mass=sum(at.get_masses())
		mass=mass*amu_to_g
		#print("mass: ",mass)
		
		volume=at.get_volume()
		volume=volume*ang3_to_cm3
		#print("volume: ",volume)
		
		density=mass/volume
		
		#print("density: ",density)
		
		#break
		
		data3['density'].append(density)
		data3['energy'].append(eng)
		
	pandas.DataFrame.from_dict(data3).to_csv("density_energy_data.csv")	
	fig2, ax1=plt.subplots()
	colour='k'
	ax1.set_xlabel("density (gcm^-1)")
	ax1.set_ylabel("Energy (eV/atom)")
	ax1.scatter(data3['density'],data3['energy'],color=colour,s=0.75)
	fig2.tight_layout()
	plt.savefig("density_energy_plot.png",dpi=600)
	
	del data, fig2, ax1
	
#search=2
#plot_graph()
