import os
import sys

#elements=['Sr','Ti','Y','O']

def assemble_spp(elements,output='lib.lib',spp_path=r'C:\Users\cc0u5\Documents\SPP_library\SPP'):
	#os.environ['SPP_PATH']=spp_path
	#print(os.environ['SPP_PATH'])
	#First enumerate the pairs of potentials:
	pairs=[]
	for i in elements:
		for j in elements:
			pair=[i,j]
			if not pair in pairs:
				pairs.append(pair)
				
	o=open(output,'w')
	
	for i in pairs:
		try:
			f=open(spp_path+"/"+i[0].upper()+"-"+i[1].upper()+"/"+i[0].upper()+"-"+i[1].upper()+".POT",'r').readlines()
		except:
			f=open(spp_path+"/"+i[1].upper()+"-"+i[0].upper()+"/"+i[1].upper()+"-"+i[0].upper()+".POT",'r').readlines()
		o.write("\n")
		for j in f:
			o.write(j)
			
	o.close()
	
	
#assemble_spp(elements)