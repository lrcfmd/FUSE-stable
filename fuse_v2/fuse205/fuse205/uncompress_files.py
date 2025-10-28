import os,sys,pickle,zipfile,shutil

#function to unpack zip archives from a previously compressed / completed FUSE run

def unzip_archive(zip_path, extract_to):
	os.makedirs(extract_to, exist_ok=True)
	with zipfile.ZipFile(zip_path, 'r') as zipf:
		zipf.extractall(extract_to)

def uncompress():

	if os.path.isfile('structures.zip'):
		unzip_archive('structures.zip', 'structures')
	
	if os.path.isfile('restart.zip'):
		unzip_archive('restart.zip', 'restart')
	
	shutil.copytree("restart","backup")
	
	if os.path.isfile('reference_structures.zip'):
		unzip_archive('reference_structures.zip', 'reference_structures')
	
	if os.path.isfile('gnboss.zip'):
		unzip_archive('gnboss.zip', 'gn-boss')
	
	