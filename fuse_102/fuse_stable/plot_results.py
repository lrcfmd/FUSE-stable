import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas
import sys
import numpy as np

def plot_graph(search=1):
	try:
		data=pandas.read_csv("graph_output.csv")
		#print (data)
		n=0
		for i in range(len(data['step'])):
			if data['current_energy'][i] == 0:
				n+=1
		
		fig1, ax1=plt.subplots()
		colour='tab:red'
		ax1.set_xlabel('Structures')
		ax1.set_ylabel('Energy (eV/atom)',color='tab:red')
		ax1.scatter(data['step'],data['energies'],color=colour,s=0.75)
		ax1.plot(data['step'][n:],data['current_energy'][n:],color='tab:green',linewidth=1)
		ax1.set_ylim([ min(data['energies'])-0.025, min(data['energies'])-min(data['energies'])*0.02 ])
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
			ax2.plot(data['step'],data['temp'],color='k',linewidth=0.75)
			ax2.tick_params(axis='y',labelcolor='k')
		fig1.tight_layout()
		#plt.show()
		plt.savefig('run_output.png',dpi=600)
	except:
		pass
#search=2
#plot_graph()
