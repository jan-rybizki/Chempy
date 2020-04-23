"""
Some utility routines copied from Alexander Heger's utils.py
"""

from os import SEEK_END
from functools import reduce
from numpy import iterable, ndarray
from operator import mul

def prod(seq):
    """Product of a sequence."""
    if not iterable(seq):
        seq = (seq,)
    return reduce(mul, seq, 1)

def cumsum0_enum_range_iter(list,seed=0):
    """Iterator for cumulative sum and counter triples for slicing starting at seed (default 0)."""
    t1 = seed
    for i,l in enumerate(list):
        t0 = t1
        t1 += l
        yield i,t0,t1

def xz_file_size(filename):
    """
    Return file size of xz xompressed files.

    The file format is documented at
    `http://tukaani.org/xz/xz-file-format.txt`
    """
    def decode(buffer, index = 0):
        """
        Decode variable length integers from buffer.
        """
        i = 0
        size = 0
        b = 0x80
        while b & 0x80 > 0:
            # b, = struct.unpack('B',buffer[index+i])
            b = buffer[index+i]
            size += (b & 0x7f) << (7 * i)
            i += 1
        return size, index+i

    with open(filename, 'rb') as f:
        f.seek(-8, SEEK_END)
        bkwd_size = ndarray((), dtype = "<u4", buffer = f.read(4))
        # print(bkwd_size)

        # 12 byte footer (4-CRC,4-size,2-flags,2-YZ)
        # 4 * (backward_size + 1) is start of index
        # index starts with 0x00 flag, last 4 byte arce CRC
        f.seek(-12 - 4 * (bkwd_size + 1), SEEK_END)
        buffer = f.read(4 * bkwd_size)
        index = 1
        num, index = decode(buffer, index)

        file_size = 0
        for i in range(num):
            # read pairs of upad_size, ucomp_size
            size, index = decode(buffer, index)
            size, index = decode(buffer, index)
            file_size += size
        return file_size
