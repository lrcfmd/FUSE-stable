"""
This module contains functionality for reading and writing an ASE
Atoms object in GULP format.

"""

from ase.geometry import *
import os

def read_gulp(filename='ase-gulp.gin'):
    """Import GULP type file.

    Reads unitcell, atom positions, charges, atom types
    Note: Constraints not currently supported
    """
 
    from ase import Atoms, Atom
    from ase.constraints import FixAtoms, FixScaled
    from ase.data import chemical_symbols
    from ase.spacegroup import crystal
    import numpy as np

    if isinstance(filename, str):
        f = open(filename)
    else: # Assume it's a file-like object
        f = filename

    #Read file
    data = f.readlines()
    #Find spacegroup, if not present default is P 1
    space_group = '1'
    for n,line in enumerate(data):
        if line.lstrip()[0:5]=='space':
            if len(line.split())>1:
                space_group = line.split()[1:].join().strip()
            else:
                space_group = data[n+1].strip()
    if space_group.isdigit():
        space_group = int(space_group)
    #Find unit cell, could be as 3x3 matrix or as parameters
    matrix = True
    cell = []
    for n,line in enumerate(data):
        if line.lstrip()[0:4]=='vect':
            for i in range(3):
                cell.append([float(data[n+i+1].split()[0]),
                             float(data[n+i+1].split()[1]),
                             float(data[n+i+1].split()[2])])
            cell = np.array(cell)
            if len(line.split())>1 and line.split()[1]=='au':
                cell = cell*ase.units.Bohr
        if line.lstrip()[0:4]=='cell':
            matrix = False
            if len(line.split())>2:
                if line.split()[1].isalpha():
                    for i in range(2,8):
                        cell.append(float(line.split()[i]))
                    cell = np.array(cell)
                    if line.split()[1]=='au':
                        cell = cell*ase.units.Bohr
                else:
                    for i in range(1,7):
                        cell.append(float(line.split()[i]))
                    cell = np.array(cell)
            else:
                for i in range(6):
                    cell.append(float(data[n+1].split()[i]))
                cell = np.array(cell)
                if len(line.split())>1 and line.split()[1]=='au':
                    cell = cell*ase.units.Bohr
    #Find atom positions, could be as fractional or cartesian coordinates
    coords = []
    symbols = []
    charges = []
    for n,line in enumerate(data):
        if line.lstrip()[0:4]=='frac' or line.lstrip()[0:4]=='cart':
            cartesian = line.find('cart')!=-1
            bohr_radius = line.find(' au')!=-1
            iatom = 0
            tmp_line = data[n+iatom+1]
            length_of_line=len(tmp_line.split())
            while len(tmp_line.split())==length_of_line:
                tmp_data = tmp_line.split()
                if tmp_data[1][0]=='c':
                    if tmp_data[0].isdigit():
                        symbols.append(ase.data.chemical_symbols[int(tmp_data[0])])
                    else:
                        symbols.append(tmp_data[0])
                    coord=[]
                    for icart in range(2,5):
                        if tmp_data[icart].find('/')!=-1:
                            coord.append(float(tmp_data[icart].split('/')[0])/
                                         float(tmp_data[icart].split('/')[1]))
                        else:
                            coord.append(float(tmp_data[icart]))
                    coords.append(coord)
                    try:
                        charges.append(float(tmp_data[5]))
                    except:
                        iatom = iatom+1
                        try:
                            tmp_line = data[n+iatom+1]
                        except:
                            break
                        continue
                elif tmp_data[1][0]!='s':
                    if tmp_data[0].isdigit():
                        symbols.append(ase.data.chemical_symbols[int(tmp_data[0])])
                    else:
                        symbols.append(tmp_data[0])
                    coord=[]
                    for icart in range(2,5):
                        if tmp_data[icart].find('/')!=-1:
                            coord.append(float(tmp_data[icart].split('/')[0])/
                                         float(tmp_data[icart].split('/')[1]))
                        else:
                            coord.append(float(tmp_data[icart]))
                    coords.append(coord)
                    try:
                        charges.append(float(tmp_data[4]))
                    except:
                        iatom = iatom+1
                        try:
                            tmp_line = data[n+iatom+1]
                        except:
                            break
                        continue
                iatom = iatom+1
                try:
                    tmp_line = data[n+iatom+1]
                except:
                    break
    coords = np.array(coords)
    #If coordinates were given as cartesian convert to fractional/scaled
    if not matrix:
        cell = cellpar_to_cell(cell)
    if cartesian:
        if bohr_radius:
            coords = coords*ase.units*Bohr
        coords = np.linalg.solve(cell.T, coords.T).T
        for i in range(3):
            coords[:, i] %= 1.0
            coords[:, i] %= 1.0
    charges=np.array(charges)
    #Check that symbols are just atomic symbols without numbers
    for isymbol in range(len(symbols)):
        if not symbols[isymbol].isalpha():
            for character in symbols[isymbol]:
                if not character.isalpha():
                    symbols[isymbol] = symbols[isymbol].strip(character)
    atoms=crystal(symbols=symbols,basis=coords,spacegroup=space_group,cell=cell)
    f.close()
    return atoms

def read_gulp_out(filename='ase-gulp.gout'):
    """Import gulp output file.

    Reads unitcell, atom positions, energies, and forces from the output file.
    Reads in final positions in an atoms object

    - at the moment the function works for GULP files made after running with
      ase, that is the initial structure is defined in the full unit cell and 
      there is only one structure in the input file
    - does not (yet) read any constraints
    """
    import os
    import numpy as np
    from ase.calculators.singlepoint import SinglePointCalculator
    from ase import Atoms, Atom

    if isinstance(filename, str):
        f = open(filename)
    else: # Assume it's a file-like object
        f = filename
    data = f.readlines()
    ncores_shells  = 0
    atoms = Atoms(pbc = True)
    energy = 0
    forces = None
    final = False
    #First read in final state if present
    for n,line in enumerate(data):
        if 'Total number atoms/shells' in line:
            ncores_shells=int(data[n].split()[4])
        if 'Final Cartesian lattice vectors (Angstroms) :' in line:
            cell = []
            for i in range(3):
                temp = data[n+2+i].split()
                cell += [[float(temp[0]), float(temp[1]), float(temp[2])]]
            atoms.set_cell(cell)
        if 'Final fractional coordinates of atoms :' in line:
            coords=[]
            for iatom in range(ncores_shells):
                temp    = data[n+6+iatom].split()
                #Only include cores
                if temp[2]=='c':
                    atoms  += Atom(temp[1],[0.,0.,0.])
                    coords.append([float(temp[3]),float(temp[4]),float(temp[5])])
            natoms=len(atoms)
            final=True
        if 'Final internal derivatives :' in line:
            forces = []
            for iatom in range(ncores_shells):
                temp    = data[n+6+iatom].split()
                #Only include cores
                if temp[2]=='c':
                        try:
                            forces += [[float(temp[3]),float(temp[4]),float(temp[5])]]
                        except (TypeError,ValueError):
                            forces = None
        if 'Final energy' in line:
            try:
                energy = float(data[n].split()[3])
            except ValueError:
                energy = data[n].split()[3]
    #If we didn't find the final information in that format look for other data
    if not final:
        for n,line in enumerate(data):
            if 'Cartesian lattice vectors (Angstroms) :' in line:
                cell = []
                for i in range(3):
                    temp = data[n+2+i].split()
                    cell += [[float(temp[0]), float(temp[1]), float(temp[2])]]
                atoms.set_cell(cell)
                final=True
            if 'Fractional coordinates of asymmetric unit :' in line:
                coords=[]
                for iatom in range(ncores_shells):
                    temp    = data[n+6+iatom].split()
                    #Only include cores
                    if temp[2]=='c':
                        atoms  += Atom(temp[1],[0.,0.,0.])
                        coords.append([float(data[n+6+iatom][18:27]),
                                       float(data[n+6+iatom][30:39]),
                                       float(data[n+6+iatom][42:51])]) 
                natoms=len(atoms)
            if 'Total lattice energy' in line and 'eV' in line:
                try:
                        energy = float(data[n].split()[4])
                except ValueError:
                        energy = data[n].split()[4]
    atoms.set_scaled_positions(coords)
    try:
	    atoms.set_calculator(SinglePointCalculator(energy,forces,None,None,atoms))
    except(AttributeError):
            atoms.set_calculator(None)
    return atoms

def write_gulp(filename, images, keywords=['sing',], options=None, shells = None, symbols = None):
    """Method to write GULP input files.
           keywords: Keywords with no parameters which appear on first line of 
                     GULP input file (Default is sing)
           options: Options with parameters. This should be a list of options 
                    with all their parameters
           shells: A dictionary object with atom labels as keys and charges as 
                   values for each atomic species to have a shell attached
           symbols: A list of symbols for each atom if they are not wanted to be
                    the same as in the original atoms object

    """
    
    import numpy as np
    from ase.constraints import FixAtoms, FixScaled

    if isinstance(filename, str):
        f = open(filename, 'w')
    else: # Assume it's a 'file-like object'
        f = filename
    
    if isinstance(images, (list, tuple)):
        if len(images) > 1 and symbols:
            raise RuntimeError("Don't know how to save more than "+
                               "one image to GULP input with specific symbols")
    else:
        images=[images,]

    #Write out keywords line
    for keyword in keywords:
        f.write("%s "%keyword)    
    f.write("\n\n")

    #Write out information for each atoms object in images
    for atoms in images:
        f.write('# number of atoms = %d\n\n' % len(atoms)) 
        # writing unitcell vectors 
        f.write('vectors\n')                              
        for vec in atoms.get_cell():                            
            f.write('  %11.6f %11.6f %11.6f\n'% tuple(vec))
        f.write("\n")
        # writing core atomic co-ordinates 
        f.write('fractional\n')
        if not symbols:
                symbols = atoms.get_chemical_symbols()
        try:
                charges=atoms.get_initial_charges()
        except(NotImplementedError,RuntimeError):
                charges=None
        if charges==None:
            charges=[]
            for atom in atoms:
                charges.append(0.)
        positions=atoms.get_scaled_positions()
        for s,(x, y, z),c in zip(symbols, atoms.get_scaled_positions(), charges):
            if shells != None:
                if shells.has_key(s):
                    c=c-shells[s]
            format_text='%s core %10.6f %10.6f %10.6f %10.6f\n'
            params=(s, x, y, z, c)
            f.write(format_text % params)
        ### writing shell atoms ###
        if shells !=None:
            for s,(x, y, z) in zip(symbols, atoms.get_scaled_positions()):
                if shells.has_key(s):
                    format_text='%s shell %10.6f %10.6f %10.6f %10.6f\n'
                    params=(s, x, y, z, shells[s])
                    f.write(format_text % params)
        f.write("\n")

    #Write out options. For each option write name followed by parameters on a new line
    if options:
        for option in options:
            f.write("%s\n"%option[0])
            for parameter in option[1:]:
                f.write("%s "%str(parameter))
            f.write("\n\n")

    if type(filename) == str:
        f.close()
