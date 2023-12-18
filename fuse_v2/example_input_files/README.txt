The four example input files here are provided to repeat the calculation setups presented in the paper (quoted from the calculation setup in the paper):

1. ”fgen”: As our baseline experiment, FUSE2 is run without using the ML model outlined in 2.1, all crystal
structure are therefore generated as in the original version of the code, using random the original unit cell
selection and sub-module motifs to populate the unit cell. Local optimisation is then only performed using
VASP. This experiment is similar to structure prediction runs using the original implementation of FUSE.

2. ”mlgen”: The initial population of crystal structures is generated using the ML model in 2.1. The ML
generated structures are ranked using SPPs and the top x structures selected and broken down into their
sub-modules. The remaining ML structures which do not form the initial population of structures then form
a pool of structures which may be introduced into the search via moves 13. and 14. outlined above. Local
optimisation is then only performed using VASP.

3. ”fgen-SPP”: As experiment 1, but local optimisation is performed in two stages: 1) structures are locally
optimised using SPPs and 2) the SPP optimised structure is then re-optimised using VASP, with the final
energy taken from VASP.

4. ”mlgen-SPP”: As experiment 2, but local optimisation is performed in two stages: 1) structures are locally
optimised using SPPs and 2) the SPP optimised structure is then re-optimised using VASP, with the final
energy taken from VASP.

Once you have chosen the calculation you wish to perform, copy the appropriate example input file & edit the required input parameters.

All of the example inputs have been configured for starting a fresh calculation for Ca3Ti2O7.