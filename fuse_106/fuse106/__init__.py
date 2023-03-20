"""Flexible Unit Structure Engine (FUSE)"""

from distutils.version import LooseVersion

import ase
import numpy



from fuse106.assemble_structure import *
#from fuse106.bond_table import *
from fuse106.create_random_instructions import *
from fuse106.create_random_string import *
from fuse106.error_check_structure import *
from fuse106.get_distances import *
from fuse106.gulp import *
from fuse106.make_basin_move import *
from fuse106.plot_results import *
from fuse106.possible_solutions import *
from fuse106.run_fuse import *
from fuse106.run_gulp import *
from fuse106.run_vasp import *
from fuse106.generate_random_structure import *
#from fuse.write_structures import *

__all__ = ['get_distances','assemble_structure','assemble_structure','return_factors','create_random_instructions',
'create_random_string','error_check_structure','read_gulp','read_gulp_out','write_gulp',
'make_basin_move','make_ga_move','plot_graph','cube_function','tetragonal_function','orthorhombic_function',
'possible_solutions','run_fuse','run_gulp','cellpar','run_vasp']

__version__ = '1.06'


if LooseVersion(np.__version__) < '1.9':
    # Make isinstance(x, numbers.Integral) work also for np.intxx:
    import numbers
    numbers.Integral.register(np.integer)
    del numbers
