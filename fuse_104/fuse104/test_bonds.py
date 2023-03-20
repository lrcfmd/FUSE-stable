from ase import *
import numpy
import sys
from ase.visualize import *

#atoms = read("0.cif")
#cations = ['Y','Ti','Ba']
#anions = ['O']
#charges = [3,4,2]
#ap = 4.95

#lib={'Y':{3:[6,9]},'Ti':{4:[4,8]},'Ba':{2:[6,12]}} # start of information from shannon database, format = 1st key: atomic symbol, 2nd key(s): charge states, then list containing min bonds then max bonds

def test_bonds(atoms='',cations='',anions='',charges='',ap='',lib='',system_type='',count_bonds=False): 
	if system_type=='ionic':
		anion_rad=[]
		#print(anions)
		for i in range(len(anions)):
			temp=lib[anions[i]]
			keys=list(temp.keys())
			if len(temp) == 1:
				anion_rad.append(temp[keys[0]][-1])
			if len(temp) > 1:
				for j in range(len(temp)):
					anion_rad.append(temp[keys[j]][-1])
		anion_rad=max(anion_rad)
		#count_bonds=True
		if count_bonds==True:
			print("########## Starting structure ##########")
			print("number of atoms: ",str(len(atoms)))
		ci=[]
		ai=[]
		#atoms.rattle(0.025)
		for i in range(len(atoms)):
			if atoms[i].symbol in cations:
				ci.append(i)
		
		atoms=atoms.repeat([2,2,2])
		for i in range(len(atoms)):
			if atoms[i].symbol in anions:
				ai.append(i)
		wrong=0
		for i in range(len(ci)):
			cat=atoms[ci[i]]
			try:
				data=lib[atoms[ci[i]].symbol][charges[cations.index(atoms[ci[i]].symbol)]]
				br=list(range(data[0]-1,data[1]+2))
			except KeyError:
				min_bonds=[]
				max_bonds=[]
				radius=[]
				temp=[]
				for j in  range(len(list(lib[atoms[ci[i]].symbol].keys()))):
					temp=lib[atoms[ci[i]].symbol][list(lib[atoms[ci[i]].symbol].keys())[j]]
					min_bonds.append(temp[0])
					max_bonds.append(temp[1])
					radius.append(temp[2])
				data=[min(min_bonds),max(max_bonds),sum(radius)/len(radius)]
				br=list(range(data[0]-1,data[1]+2))
			bond_dist=(data[-1]+anion_rad)*1.25
			if count_bonds==True:
				print("bond_distance ", str(bond_dist))
			nbonds=0
			temp1=list(atoms.get_distances(ci[i],ai,mic=True))
			#print(temp1)
			try:
				temp1.remove(0.)
			except:
				pass
			###### to restore the bond test used in v. 1.02 chage the distance comparison to: if temp1[j] <= ap*0.785
			for j in range(len(temp1)):
				if temp1[j] <= bond_dist:
					nbonds+=1
			#print(temp1)
			#for j in range(len(ai)):
			#	ani=atoms[ai[j]]	
			#	bond=atoms.get_distance(ci[i],ai[j],mic=True)
			#	if bond <= ap*0.72:
			#		nbonds += 1
			#		#if atoms[j].symbol in cations:
			#		#	wrong+=100
			#print(cations)
			#print(charges)
			#sys.exit()
			#print(str(atoms[ci[i]].symbol),end=" ")
			#print(nbonds)
			try:
				data=lib[atoms[ci[i]].symbol][charges[cations.index(atoms[ci[i]].symbol)]]
				br=list(range(data[0]-1,data[1]+2))
			except KeyError:
				min_bonds=[]
				max_bonds=[]
				radius=[]
				temp=[]
				for j in  range(len(list(lib[atoms[ci[i]].symbol].keys()))):
					temp=lib[atoms[ci[i]].symbol][list(lib[atoms[ci[i]].symbol].keys())[j]]
					min_bonds.append(temp[0])
					max_bonds.append(temp[1])
					radius.append(temp[2])
				data=[min(min_bonds),max(max_bonds),sum(radius)/len(radius)]
				br=list(range(data[0]-1,data[1]+2))
			#print(data)
			if count_bonds==True:
				print (cat.symbol, " ",nbonds)
			if not nbonds in br:
				wrong += 1
		#print ("total errors", wrong)
		if count_bonds==True:
			write("test_structure.cif",atoms)
		#sys.exit()
		try:
			return float(wrong)/float(len(ci))
		except:
			return 1.
		
	if system_type=='neutral':
		ci=[]
		ai=[]
		#atoms.rattle(0.025)
		for i in range(len(atoms)):
			if atoms[i].symbol in cations:
				ci.append(i)
			if atoms[i].symbol in anions:
				ai.append(i)
		wrong=0
		for i in range(len(ci)):
			cat=atoms[ci[i]]
			nbonds=0
			for j in range(len(ai)):
				ani=atoms[ai[j]]	
				bond=atoms.get_distance(ci[i],ai[j],mic=True)
				if bond <= ap*0.785:
					nbonds += 1
			
			min_bonds=[]
			max_bonds=[]
			radius=[]
			temp=[]
			for j in  range(len(list(lib[atoms[ci[i]].symbol].keys()))):
				temp=lib[atoms[ci[i]].symbol][list(lib[atoms[ci[i]].symbol].keys())[j]]
				min_bonds.append(temp[0])
				max_bonds.append(temp[1])
				radius.append(temp[2])
				
			data=[min(min_bonds),max(max_bonds),sum(radius)/len(radius)]
			br=list(range(data[0]-1,data[1]+2))
			#print cat.symbol, " ",nbonds
			if not nbonds in br:
				wrong += 1
		#print ("total errors", wrong)
		#write("test_structure.cif",atoms)
		#sys.exit()
		try:
			return float(wrong)/float(len(ci))
		except:
			return 1.

#wrong_fract=test_bonds(atoms=atoms,cations=cations,anions=anions,charges=charges,ap=ap,lib=lib)
#print wrong_fract


