#set the string in os.system to be the required run command for vasp
import os
exitcode = os.system('timeout --signal=9 60m mpirun -np 36 vasp_std >vasp_output.txt')

