#! /bin/env python3

"""
A collection of tool using fortran IO library
"""

try:
    from .reader import FortranReader
    from .writer import FortranWriter
except:
    from fortranfile.reader import FortranReader
    from fortranfile.writer import FortranWriter

# TODO - accept lists, slices, and append to file
def extract_record(filein, fileout, nr = None):
    """extract record (0-based)"""
    with FortranReader(filein) as f:
        if nr is not None:
            f.skip(nr)
            f.load()
        else:
            # load last record
            while not f.eof():
                f.load()
        data = f.get_data()
        print('Extracting record {} ({} bytes)'.format(f.rpos-1, len(data)))
        byteorder = f.byteorder
    with FortranWriter(fileout, byteorder = byteorder) as f:
        f.write_data(data)

import sys

if __name__ == "__main__":
    extract_record(*sys.argv[1:])
