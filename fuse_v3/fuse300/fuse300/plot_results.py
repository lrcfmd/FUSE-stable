import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas
import sys
import numpy as np
import glob
import os
from ase.io import *
from ase import Atoms

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
	plt.close(fig1)
	del fig1,ax1
	try:
		del ax2
	except:
		pass
	
	#create plot of used moves distribution
	used_moves={1:0,2:0,3:0,4:0,5:0,6:0,7:0,8:0,9:0,10:0,11:0,12:0,13:0,14:0,15:0}
	
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
	
	plt.close(fig2)
	del data, fig2, ax1
	
	# build histogram of energy distrubution relative to global minimum
	bin_size = 5 # meV/atom
	raw=data3['energy']
	#print(raw)
	raw.sort()
	#print("2")
	raw=np.array(raw,dtype='float')
	#print("1")
	rel=raw-min(raw)
	rel=rel*1000
	#print("hello")
	#print(rel)
	upper_limit=1000
	max_bin=0
	
	for entry in rel:
		while entry > max_bin:
			max_bin+=bin_size
	
	max_bin+=bin_size+bin_size
	
	if max_bin > upper_limit:
		max_bin = upper_limit
	
	plot_data={'x':[],'y':[]}
	bin=0
	while bin < max_bin:
		plot_data['x'].append(bin)
		plot_data['y'].append(0)
		bin+=bin_size
	
	bins=np.arange(0-bin_size,max_bin,bin_size)

	indicies=np.digitize(rel,bins)-1
	
	counts=np.bincount(indicies,minlength=len(bins))
	
	histogram={'bin':bins.tolist(),'counts':counts.tolist()}
	pandas.DataFrame.from_dict(histogram).to_csv("histogram_data.csv",index=False)
	
	# plot it
	fig3,ax1=plt.subplots()
	ax1.set_xlabel("Relative energy (meV/atom), Start of bin")
	ax1.set_ylabel("Counts (a.u.)")
	#ax1.scatter(histogram['bin'],histogram['counts'],s=0.75)
	ax1.plot(histogram['bin'],histogram['counts'],'-o', color='black',markersize=0.5,linewidth=0.5)
	ax1.fill_between(histogram['bin'], histogram['counts'], color='green', alpha=0.4)
	
	fig3.tight_layout()
	plt.savefig("energy_distribution",dpi=600)
	plt.close(fig3)
	del fig3,ax1	
	
	# create a bar chart for the number of atoms used in the calculation
	#print(type(initial_population))
	#print(initial_population)
	data={}
	
	for i in list(initial_population.keys()):
		at=initial_population[i]['atoms'].copy()
		if not len(at) in list(data.keys()):
			data[len(at)] = [1]
		else:
			data[len(at)][0] +=1
	
	plot_data={'x':[],'y':[]}
	for i in list(data.keys()):
		plot_data['x'].append(i)
		plot_data['y'].append(data[i][0])
	
	histogram={'bin':plot_data['x'],'counts':plot_data['y']}
	pandas.DataFrame.from_dict(histogram).to_csv("n_atoms_histogram_data.csv",index=False)
	#print(histogram)

	fig4,ax4=plt.subplots()
	ax4.bar(histogram['bin'],histogram['counts'],edgecolor='black',linewidth=0.75,
		alpha=0.9,color='gold',width=0.9)
	
	ax4.set_xlabel('N atoms')
	ax4.set_ylabel('counts')
	fig4.tight_layout()
	fig4.savefig("atom_distribution.png",dpi=600)
	plt.close(fig4)
	del fig4,ax4
	#plot energy vs. number of atoms
	data={}
	
	plot_data={'x':[],'y':[]}
	for i in list(data.keys()):
		plot_data['x'].append(i)
		plot_data['y'].append(data[i][0])
		
	for i in list(initial_population.keys()):
		at=initial_population[i]['atoms'].copy()
		if initial_population[i]['converged']==True:
			plot_data['x'].append(len(at))
			plot_data['y'].append(initial_population[i]['energy'])	
	
	histogram={'atoms':plot_data['x'],'energy':plot_data['y']}
	pandas.DataFrame.from_dict(histogram).to_csv("atoms_energy_data.csv",index=False)
	#print(histogram)

	fig5,ax5=plt.subplots()
	ax5.scatter(histogram['atoms'],histogram['energy'],
		alpha=0.9,color='black',s=0.75)
	
	ax5.set_xlabel('N atoms')
	ax5.set_ylabel('energy (eV/atom)')
	fig5.tight_layout()
	fig5.savefig("atom_v_energy.png",dpi=600)
	plt.close(fig5)
	del fig5,ax5

#search=2
#plot_graph()
