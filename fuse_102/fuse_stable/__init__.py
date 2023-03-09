"""Flexible Unit Structure Engine (FUSE)"""

from distutils.version import LooseVersion

import ase
import numpy

from fuse_stable.assemble_structure import *
from fuse_stable.create_random_instructions import *
from fuse_stable.create_random_string import *
from fuse_stable.error_check_structure import *
from fuse_stable.get_distances import *
from fuse_stable.gulp import *
from fuse_stable.make_basin_move import *
from fuse_stable.plot_results import *
from fuse_stable.possible_solutions import *
from fuse_stable.run_fuse import *
from fuse_stable.run_gulp import *
from fuse_stable.run_vasp import *

__all__ = ['get_distances','assemble_structure','assemble_structure','return_factors','create_random_instructions',
'create_random_string','error_check_structure','read_gulp','read_gulp_out','write_gulp',
'make_basin_move','plot_graph','cube_function','tetragonal_function','orthorhombic_function',
'possible_solutions','run_fuse','run_gulp','cellpar','run_vasp']

__version__ = '1.04'


if LooseVersion(np.__version__) < '1.9':
    # Make isinstance(x, numbers.Integral) work also for np.intxx:
    import numbers
    numbers.Integral.register(np.integer)
    del numbers
