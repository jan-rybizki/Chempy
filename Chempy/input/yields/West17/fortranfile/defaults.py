"""
This class sets the defaults for Fortran I/O
"""

import numpy as np

class FortranSpecs(object):
    """
    Specification of FORTRAN file characteristics
    """

    def _set_reclen(self, reclen):
        if reclen == 8:
            dtype = np.uint64
        elif reclen == 4:
            dtype = np.uint32
        else:
            raise Exception('Invald reclen.')
        self.reclen_dtype = np.dtype(dtype)
        self.fortran_reclen = self.reclen_dtype.itemsize
