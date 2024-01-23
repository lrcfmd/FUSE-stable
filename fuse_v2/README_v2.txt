Last Updated: 22nd January 2024

This is the current stable version of FUSE v.2, and is provided "as is" under the GNU public licence:


Version 2.02 and later are detailed in the following pre-print:

https://chemrxiv.org/engage/chemrxiv/article-details/658065cb9138d23161f68ac6

We have tested FUSE v.2 primarily on Linux Desktops and clusters.

We recommend setting up FUSE v.2 using Anaconda with python 3.9 or newer, as this allows you to easily create the virtual environment needed to setup the ML model used for structure generation.

********************************************************************************

FUSE v.2 is dependant on the following packages:

The atomic simulation environment (ase)
Pandas
spglib
func_timeout

all are available using pip

FUSE also uses external chemistry codes in order to perform energy calculations, which you will need acess to in order to perform CSP calculations.

Currently, GULP, VASP are supported, with Quantum Espresso implemented, but this has not been tested extensively.
In principle, any calculator which is supported by ase should be useable (https://wiki.fysik.dtu.dk/ase/).
If there is a currently un-supported calculator that you would like to use, please to get in touch!

We have additionally implemented the use of the CHGNet ML potential published here by B. Deng et al. :

https://www.nature.com/articles/s42256-023-00716-3

This can be installed using pip:

python -m pip install chgnet

If you wish to use CHGNet in your calculations, you will need to download the master version of ase from here:

https://gitlab.com/ase/ase/-/tree/master/ase?ref_type=heads

You can install this version by downloading / unpacking the zip archive & installing with pip in the root ase directory:

python -m pip install . 

If you have previously installed ase on your system you will need to append the above with "--force-reinstall"

We are currently testing how to best use CHGNet within FUSE, and will publish example input scripts for this once we have completed testing.

you can then install FUSE by typing:

python -m pip install in the root python directory.

**********************************************************************************

Installing the gnboss ML model.

The ML structure generater is replacated from the model presented in this paper by G. Cheng et al.:

https://www.nature.com/articles/s41467-022-29241-4

to install, first unpack the "gnboss_template.zip" archive and go into the first directory.

create a fresh anaconda environment:

conda create --name gnboss python=3.6.13

activate the environment:

conda activate gnboss

then install the required packages:

python -m pip install -r requirements.txt

then deactivate the model by typing:

conda deactivate.

**********************************************************************************

To use the statistical proxy potential library provided, unzip the "SPP_light.zip" and take note of path to the unpacked directory.

The statistical proxy potentials are implemented for use with GULP, and are published in this paper by D. Antypov et al.:

https://chemrxiv.org/engage/chemrxiv/article-details/65292c4a8bab5d20554598dd

**********************************************************************************

Finally, to use all of the features for FUSE v.2, it is reccomended to set the following environment variables:

for GULP:
"ASE_GULP_COMMAND" : the path where your gulp executable can be found, followed by "gulp < gulp.gin > gulp.got"
"GULP_LIB": This can be set to just empty quotation marks.

for VASP:
"VASP_PP_PATH": The path to the directory containing vasp pseudo-potentials
"VASP_SCRIPT": The path to the python script for running vasp. For the example input files provieded here, this just needs to be "vasp.py"

for using the gnboss ML model:
"CONDA": that path to the location of your conda executable, this is required to run the ML structure geneartion model in it's seperate conda virtual environment. this is normally located in the directory that you have installed anaconda in the "bin" directory.
"GNBOSS": the path to the virtual environment that you have created above for the gnboss model, this is typically located in the "envs" folder in the loaction that you have installed anaconda.
"GNBOSS_TEMP": The path to the location of the gnboss template that you have unpacked above.

for using statistical proxy potentials (this also requires GULP):
"SPP_PATH": The path to the statistical proxy potential directory unpacked above, so that FUSE can asseble any required potential sets.




