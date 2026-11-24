import os
import sys

#elements=['Sr','Ti','Y','O']

def assemble_spp(elements,output='lib.lib',spp_path=r'C:\Users\cc0u5\Documents\SPP_library\SPP'):
	
	pairs=[]
	
	for i in elements:
		for j in elements:
			pair=[i,j]
			if not pair in pairs:
				pairs.append(pair)
				
	o=open(output,'w')
	
	found_pairs=[]
	
	for i in pairs:
		if i in found_pairs:
			continue
		
		j=i.copy()
		j.reverse()
		
		if j in found_pairs:
			continue
		
		
		made=False
		try:
			f=open(spp_path+"/"+i[0].upper()+"-"+i[1].upper()+"/"+i[0].upper()+"-"+i[1].upper()+".POT",'r').readlines()
			
		except:
			f=open(spp_path+"/"+i[1].upper()+"-"+i[0].upper()+"/"+i[1].upper()+"-"+i[0].upper()+".POT",'r').readlines()
			
		found_pairs.append(i)	
		
		o.write("\n")
		for x in f:
			o.write(x)
			
	o.close()
	print(found_pairs)
	
	
#assemble_spp(elements)