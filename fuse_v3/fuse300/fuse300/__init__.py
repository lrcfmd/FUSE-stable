"""Flexible Unit Structure Engine (FUSE)"""

from distutils.version import LooseVersion

import ase
import numpy

from fuse300.bond_table import *
from fuse300.extract_modules import *

__all__ = ['bond_table','extract_modules']

__version__ = '3.00'


if LooseVersion(numpy.__version__) < '1.9':
    # Make isinstance(x, numbers.Integral) work also for np.intxx:
    import numbers
    numbers.Integral.register(np.integer)
    del numbers
