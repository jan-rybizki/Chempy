"""
Classes for reading and writing UNIX unformatted FORTRAN files.
"""

# TODO
# * for general reading, load needs to be able to specify record size?

from .reader import FortranReader
from .writer import FortranWriter
