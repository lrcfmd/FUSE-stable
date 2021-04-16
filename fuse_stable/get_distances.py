from ase import *

def get_distances(new_atoms=''):
	distances=[]
	for i in range(len(new_atoms)):
		for x in range(len(new_atoms)):
			if i != x:
				distances.append(new_atoms.get_distance(i,x,mic=True))
	return distances