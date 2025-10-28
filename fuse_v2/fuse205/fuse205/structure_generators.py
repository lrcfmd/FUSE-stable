import os
import shutil
import glob
import ase
from ase.io import *
import pandas

from fuse205.run_castep import *
from fuse205.run_gulp import *
from fuse205.run_vasp import *
from fuse205.run_qe import *
from fuse205.assemble_spp import *

homedir=os.getcwd()

def run_generators(o='',generator='',generate_gn_boss_structures='',clear_previous_structures='',gn_boss_command='',gn_search='',gn_max_step='',
                   gn_template_path='',composition='',gn_zn_range='',path_to_structures='',max_atoms='',generate_airss_structures='',
                   airss_form_units='',airss_num_structures='',repose_seed='',repose_command='',repose_devmax=''):

    if not os.path.isdir(path_to_structures):
        os.mkdir(path_to_structures)
    if clear_previous_structures == True:
        shutil.rmtree(path_to_structures)
        os.mkdir(path_to_structures)

    for x in generator:
        if x == 'gnboss':
            if generate_gn_boss_structures == True:
                print("Generating structure pool using Gn-Boss ML model")
                o.write("\nGenerating structure pool using Gn-Boss ML model\n")

                if clear_previous_structures == True:
                    if os.path.isdir("gn-boss"):
                        shutil.rmtree("gn-boss")  # clear out the previous run.
                        os.mkdir("gn-boss")
                    else:
                        os.mkdir("gn-boss")
                    os.chdir("gn-boss")

                if clear_previous_structures != True:
                    if not os.path.isdir("gn-boss"):
                        os.mkdir("gn-boss")

                    os.chdir("gn-boss")

                # go fetch the template file
                # try:
                shutil.copytree(gn_template_path, '.', dirs_exist_ok=True)
                # os.system("cp -r "+gn_template_path+"\* .")
                # except:
                # os.system("cp -r "+gn_template_path+"* .")

                os.chdir("chemical_compositions")
                template = open("template.in", 'r').readlines()
                run_files = []
                # for each value of Z build the input file
                for z in range(gn_zn_range[0], gn_zn_range[-1] + 1):
                    run_file = template.copy()
                    form = 'compound = '
                    total = 0
                    for w in list(composition.keys()):
                        form += w
                        form += str(composition[w] * (z))
                        total += composition[w] * (z)
                        form += " "
                    form += "\n"
                    if total <= max_atoms:
                        run_file[2] = form

                    run_file[35] = "algorithm = " + gn_search + "\n"
                    run_file[39] = "max_step = " + str(gn_max_step) + "\n"

                    run_file2 = open("gnoa-input_" + str(z) + ".in", 'w')
                    for w in run_file:
                        run_file2.write(w)

                    run_file2.close()

                    run_files.append("gnoa-input_" + str(z) + ".in")

                os.remove("template.in")

                os.chdir("../")
                #print(os.getcwd())
                # sys.exit()
                # now go and run the calculations
                os.system(gn_boss_command)
                # collate the results in the reference structures folder
                os.chdir("results")
                r_files = glob.glob("*")
                for w in r_files:
                    if os.path.isdir(w):
                        if not w == "best_structures":
                            os.chdir(w)
                            # print(os.getcwd())
                            # sys.exit()
                            os.chdir("structures")
                            for i in glob.glob('*'):
                                os.rename(str(i),'gnboss_'+str(i))

                            # print(os.getcwd())
                            # sys.exit()
                            to_copy = glob.glob("*.cif")
                            for v in to_copy:
                                if v != 'temp.cif':
                                    shutil.copy(v, "../../../../" + path_to_structures + "/.")
                            os.chdir("../../")

                os.chdir("../../")
                os.chdir(homedir)

            # print(os.getcwd())

        if x == 'airss':
            if generate_airss_structures == True:
                print("Generating structure pool using AIRSS")
                o.write("\nGenerating structure pool using AIRSS\n")

                if clear_previous_structures == True:
                    if os.path.isdir("airss"):
                        shutil.rmtree("airss")  # clear out the previous run.
                        os.mkdir("airss")
                    else:
                        os.mkdir("airss")
                    os.chdir("airss")

                if clear_previous_structures != True:
                    if not os.path.isdir("airss"):
                        os.mkdir("airss")
                    os.chdir("airss")

                # write seed .cell file
                for i in airss_form_units:
                    species=''
                    id=''
                    for k in composition:
                        species += (str(k)+'%NUM='+str(composition[k])+',')
                        id += (str(k)+str(int(composition[k])*i))

                    with open('template_seed.cell','w') as file:
                        file.write('#SPECIES='+str(species)[:-1]+'\n')
                        file.write('#NFORM='+str(i)+'\n')
                        file.write('#MINSEP=1.0'+'\n')

                    for k in range(airss_num_structures):
                        os.system('buildcell < template_seed.cell > out.cell')
                        struc=read('out.cell')
                        write('out.cif',struc,format='cif')
                        name='airss_'+str(id)+'_'+str(k)+'.cif'
                        os.rename('out.cif', name)

                # move them to reference_structures
                os.chdir('../')
                os.chdir(path_to_structures)
                os.system('cp ../airss/*.cif .')
                os.chdir('../')

# if generator == diffcsp
# if generator etc.

def rank_gen_structures(o='',ranking='',path_to_structures='',assemble_spp_='',vasp_opts='',kcut='',produce_steps='',shel='',spp_path='',
                        use_spglib='',r_kwds='',r_gulp_opts='',r_lib='',r_calcs='',dist_cutoff='',qe_opts='',rank_structures='',
			 			gulp_command='',gulp_timeout='',n_opts='',rel='',relaxer_opts='',opt_class='',opt_device='',mode='',repose_seed='',repose_command='',
			 			repose_devmax=''):

    try:
        os.mkdir(path_to_structures)
        os.chdir(path_to_structures)

    except(FileExistsError):
        os.chdir(path_to_structures)

    # check to see that there's some files in here to use!
    r_cifs = glob.glob("*.cif")

    if len(r_cifs) != 0:

        r_results = {'file': [], 'energy': [], 'atoms': [], 'converged': []}

        if not os.path.isfile("dummy.lib"):
            fr = open("dummy.lib", 'w')
            fr.close()

        # get the required spp library
        temp = read(r_cifs[0])
        r_elements = []
        for w in range(len(temp)):
            if not temp[w].symbol in r_elements:
                r_elements.append(temp[w].symbol)

        # sys.exit()
        if assemble_spp_ == True:
            assemble_spp(r_elements, spp_path=spp_path)

        if ranking == 'gulp':
            print("Ranking generated structures using SPPs")
            o.write("\nRanking generated structures using SPPs")
        if ranking == 'chgnet':
            print("Ranking generated structures using CHGnet")
            o.write("\nRanking generated structures using CHGnet")
        if ranking == 'orb':
            print("Ranking generated structures using Orb")
            o.write("\nRanking generated structures using Orb")
        if ranking == 'mixed':
            print("Ranking generated structures using mix of calculators")
            o.write("\nRanking generated structures using mix of calculators")

        for w in range(len(r_cifs)):
            print(str(w + 1) + " of: " + str(len(r_cifs)), end='\r')
            try:
                atoms = read(r_cifs[w])
            except:
                continue
            if ranking == 'gulp':
                try:
                    atoms, energy, converged = run_gulp(atoms=atoms, shel=shel, kwds=r_kwds, opts=r_gulp_opts,
                                                        lib=r_lib, produce_steps=False, gulp_command=gulp_command,
                                                        gulp_timeout=gulp_timeout, use_spglib=use_spglib)
                except:
                    converged = False
                    energy = 1.e20

            if ranking == 'chgnet':
                from fuse205.run_chgnet import run_chgnet
                try:
                    atoms, energy, converged = run_chgnet(atoms, n_opts=n_opts, relaxer_opts=relaxer_opts,rel=rel,
                                                          opt_class=opt_class, mode=rank_structures,
                                                          opt_device=opt_device, use_spglib=use_spglib)
                except:
                    converged = False
                    energy = 1.e20

            if ranking == 'orb':
                from fuse205.run_orb import run_orb
                try:
                    atoms,energy,converged =run_orb(atoms)
                except:
                    converged = False
                    energy = 1.e20

            if ranking == 'mixed':

                try:
                    atoms, energy, converged = run_calculators(atoms=atoms, vasp_opts=
                    vasp_opts, kcut=kcut, produce_steps=None, shel=shel, use_spglib=use_spglib,
                                                               kwds=r_kwds, gulp_opts=r_gulp_opts, lib=r_lib,
                                                               calcs=r_calcs, dist_cutoff=dist_cutoff,
                                                               qe_opts=qe_opts, rel=rel,
                                                               gulp_command=gulp_command, gulp_timeout=gulp_timeout,
                                                               n_opts=n_opts, relaxer_opts=relaxer_opts,
                                                               opt_class=opt_class,
                                                               opt_device=opt_device, mode=rank_gn_structures)

                except:
                    converged = False
                    energy = 1.e20


            if ranking == 'repose':
                try:
                    write(repose_seed+".cell",atoms)	
                    #print(repose_command+" "+repose_seed+" + > repose.out")
                    os.system(repose_command+" "+repose_seed+" > repose.out")
                    output=None
                    try:
                    	  output=open("repose.out",'r').readlines()
                    except:
                    	  pass
                    energy=None
                    dev=0
                    deviation=None
                    try:
                        for z in output:
                        	if "Enthalpy:" in z:
                        			energy=z
                        			#break
                        	if "Deviation:" in z:
                        			deviation=z
                        			#break
                     		
                        energy = float(energy.split(":")[-2])
                        if deviation != None:
                            dev = float(deviation.split(":")[-1])
                        atoms=read(repose_seed+"-out.cell")
                        converged=True
                    except:
                        converged=False
                        energy=1.e20
                      
                except:
                    converged=False
                    energy=1.e20
                
                if dev > repose_devmax:
                    converged=False
                    energy=1.e20
         		

            energy = energy / len(atoms)

            r_results['file'].append(r_cifs[w])
            r_results['energy'].append(energy)
            r_results['atoms'].append(atoms)
            r_results['converged'].append(converged)

            write(r_cifs[w], atoms)

        dat = pandas.DataFrame.from_dict(r_results).sort_values(['energy'], axis=0, ascending=True)
        dat2 = dat.to_dict(orient='list')
        # print(dat2.keys())
        table = {'file': [], 'energy': []}
        for j in range(len(dat2[list(dat2.keys())[0]])):
            table['file'].append(dat2['file'][j])
            table['energy'].append(dat2['energy'][j])

        dat3 = pandas.DataFrame.from_dict(table)
        dat3.to_csv("ranking.csv", index=None)
        os.chdir('../')

    if len(r_cifs) == 0:
        os.chdir('../')
