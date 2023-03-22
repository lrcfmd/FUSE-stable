Last Updated: 22nd March 2023

This is the current stable version of FUSE, originally published in this paper: and is provided "as is" under the GNU public licence

https://pubs.rsc.org/en/content/articlelanding/2018/fd/c8fd00045j#!divAbstract
********************************************************************************
There are now three stable versions of FUSE available:

102: The python3 implementation of the code presented in the original paper linked above

104: The main change to the code is the implementation of searches with compositions without needing to specify oxidation states,
this results in FUSE treating all atoms in the same way, and is useful in systems where the number of cations greatly outnumbers 
the number of anions

106: The main changes from 104 above, is that the random structure generation has been re-written such that FUSE now samples more
evenly from the number of available formular units, where as previously it was (un-intentionally) biased towards the smaller numbers
of formula units. With the Basin hopping moves, the search routine has been altered such that when a chosen move fails to generate a valid
structure, FUSE will attempt the same move several more times before choosing a different move. In order to offer more control over move 1 (swapping
atoms between sub-modules), this has now been split into three moves (1, 8 and 9) where 1 = swap two atoms, 8 = swap 3 -> n atoms and 9 = swap all 
atoms.

********************************************************************************

FUSE is dependant on the following packages:

The atomic simulation environment (ase)
Pandas 

both are available using pip

FUSE also uses external chemistry codes in order to perform energy calculations, which you will need acess to in order to perform CSP calculations.

Currently, GULP and VASP are supported, with Quantum Espresso and CASTEP on the to-do list!
In principle, any calculator which is supported by ase should be useable (https://wiki.fysik.dtu.dk/ase/).
If there is a currently un-supported calculator that you would like to use, please to get in touch!

you can then install FUSE by typing:

python setup.py install

**********************************************************************************

running FUSE

to run FUSE, after configuring your energy calculator (see below), type:

python < [MY INPUT FILE].py > [MY OUTPUT]

**********************************************************************************

The input file:

each input file starts with the import commands "from fuse_stable import *", "import sys" and "import time".

We this download, we have included example input files in the folder "examples".

Below are the parameters which are common to any FUSE calculation:

"composition" : this parameter sets the emperical formula for the calculation, this is in the format of a dictionary, where the keys are element symbols, and the values are a 2 element list: the number of atoms in the formula and the formal charge state,
e.g. the formula SrO would be entered as: composition = {'Sr':[1,+2],'O':[1,-2]}

"rmax" : This is the parameter which controls the convergance of the basin hopping loop, FUSE counts how many structures since the current global minimum was found, once this count reaches rmax, the calculation will stop. For probe structure calculations, we will normally set this to somewhere between 500 - 1000

"iterations" : This is the number of structures FUSE will relax for the given run

"restart" : If set to FALSE, FUSE starts a fresh calculation, if set to TRUE, FUSE will attempt to restart a previous run, note, if FUSE is unable to detect / read the restart files, it will start a fresh calculation

"initial_gen" : This is the number of structures for FUSE to use in it's initial population. Typical values for this are 5 - 50

"search_gen" : This is the number of structures FUSE will generate & relax at each step of the basin hopping routine, this is typically left at 1, but we have found it useful to increase this value for systems with a large number of atoms (e.g. > 50)

"max_atoms" : This sets the maximum number of atoms which FUSE can use in a given structure, for probe structures, we will typically set this to < 70, but we have tested it for upto 300.

"imax_atoms" : This controls the maximum number of atoms which FUSE can use in ONLY the structures in the initial population. Typically, this will be equal to max_atoms, however for systems containing greater than 50 atoms, it can often be useful to set imax_atoms to 50 (or equal to 1 formula unit, if the emperical formula contains more than 50 atoms!)

"ctype" : This tells FUSE which energy calculator to use, currently supported are "gulp" and "vasp"

Note: at the bottom of every input file, there is the call "run_fuse( .... " this calls the main FUSE function and is fed all of the option from the input file.

***********************************************************************************

GULP

We have provided an example input file and interatomic potentials for running 20 atoms of SrTiO3, the interatomic potentials are taken from the original FUSE paper: https://pubs.rsc.org/en/content/articlelanding/2018/fd/c8fd00045j#!divAbstract

Note, FUSE assumes that the GULP calculator is setup as per the instructions from ase (https://wiki.fysik.dtu.dk/ase/ase/calculators/gulp.html#module-ase.calculators.gulp)

This covers the options which are specific to the GULP example input file:

note, that with the GULP calculator in FUSE, it can run GULP multiple times on the same structure, this can often be useful to effciently relax structures.
For each attempt you wish to do for each structure, you will need additional entries in "kwds" and "gulp_opts" below. In the example file supplied, we perform GULP calculations in two stages.

"kwds" :  This is a list, with each element being a string containing the keywords to use in your GULP calculation, e.g. ("opti conp")

"gulp_opts": This is a list, with each element being string of optional keywords for GULP, note, in order to correctly create a GULP input file, each option needs seperating by a \n to print it to a new line. (e.g. "\nlibrary lib2.lib \ndump temp.res")

"lib" : This is to tell GULP the name of any interatomic potential library to use for the calculation, if it is supplied as part of the options above, this can be left as ''

'shel' : This is a list of any elements which require a shell adding to them in GULP, (e.g. ['Ba']). If no shells are being used, leave as [''].

***********************************************************************************

TODO:

Provide and example input file for using FUSE with VASP (along with update this file with the VASP options!). 

Incorporate Quantum Espresso into FUSE (currently in testing stages)

Incorporate CASTEP into FUSE

***********************************************************************************

KNOWN ISSUES

 - There seems to be an intermitant issue when using Linux / MAC machines where FUSE is unable to correctly read the output from GULP calculations, this manifests as each structure being listed as "failed" and the corrisponding energy set to 1E20 / number of atoms.

