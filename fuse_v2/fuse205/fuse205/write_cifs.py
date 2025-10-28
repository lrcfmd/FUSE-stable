import os,sys,pickle,zipfile,glob
from ase.io import *

# function to write structures from fuse run as .cifs

def unzip_archive(zip_path, extract_to):
	os.makedirs(extract_to, exist_ok=True)
	with zipfile.ZipFile(zip_path, 'r') as zipf:
		zipf.extractall(extract_to)

def write_cifs():

	if not os.path.isdir("cifs"):
		os.mkdir("cifs")
		
	if not os.path.isdir("structures"):
		if os.path.isfile("structures.zip"):
			unzip_archive('structures.zip', 'structures')
			
	if os.path.isdir("structures"):
		os.chdir("structures")
		for i in glob.glob("*.p"):
			dat=pickle.load(open(i,'rb'))
			write(f"../cifs/{i[0:-2]}.cif",dat['atoms'])
		os.chdir("../")
	

		
	