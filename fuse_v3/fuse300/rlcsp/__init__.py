from distutils.version import LooseVersion

import numpy as np

from .rlcsp import *
from .state import *

__all__ = ['rlcsp.py', 'state']

__version__ = '1.10'

if LooseVersion(np.__version__) < '1.9':
    # Make isinstance(x, numbers.Integral) work also for np.intxx:
    import numbers
    numbers.Integral.register(np.integer)
    del numbers
