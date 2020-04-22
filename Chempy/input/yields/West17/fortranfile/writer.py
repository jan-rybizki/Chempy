"""
Classes for writing UNIX unformatted FORTRAN files.
"""

# TODO
# * for general reading, load needs to be able to specify record size?
# * fortranfile needs "backspace" and "truncate" functions.

import os
import sys
import gzip
import bz2
import lzma

import numpy as np

from .utils import prod
from .errors import RecordBeginningError, WriteError
from .common import _np_types, _set_method
from .defaults import FortranSpecs

#=======================================================================
# WRITER
#=======================================================================

class DataOutputBuffer(object):
    """
    Interface for writing buffered output data
    """

    default_byte_order = ">"

    sys_is_le = sys.byteorder == "little"
    native_byteorder = "<" if sys_is_le else ">"
    initial_buf_size = 2**24
    buf_grow_factor = 2
    buf_grow_limit = 2**28
    # as of 2011, even my best solid state drive will take 0.5 s to write that much

    def __init__(self,
                 byteorder = default_byte_order,
                 **kwargs):
        self._set_byteorder(byteorder = byteorder)
        self._init()


    def _init(self):
        self.pos  = 0
        self.buffer = bytearray(self.initial_buf_size)
        self.buf_size = self.initial_buf_size

    def _set_byteorder(self,
                       byteorder = default_byte_order):
        """
        set up all data types for deserved byte order
        """
        if byteorder == "=":
            byteorder = self.native_byteorder
        self.swapbyteorder = byteorder != self.native_byteorder
        self.byteorder = byteorder

    def bor(self):
        """Return whether position is beginning of record."""
        return self.pos == 0

    def assert_bor(self):
        """
        Throw exception if current position is not beginnig of record.

        This can be used to deterime whether all previous data has been written,
        i.e., as a consistency check of previous writes.
        """
        if not self.bor():
            raise RecordBeginnigError(self.pos)

    def _extend_buf(self):
        """
        Grow write buffer as specified.
        """
        self.buf_size += min(self.buf_size, self.buf_grow_limit)
        new_buffer = bytearray(self.buf_size)
        new_buffer[0:self.pos] = self.buffer[0:self.pos]
        del self.buffer
        self.buffer = new_buffer

    def _check_buf_size(self, size, offset = None):
        if offset is not None:
            if offset <= 0:
                p = size - offset
            else:
                p = self.pos + size + offset
        else:
            p = self.pos + size
        while p > self.buf_size:
            self._extend_buf()

    def skip_bytes(self, nbytes, fill = None):
        """Skip a number of empty bytes, optionally initializing with `fill`."""
        self._check_buf_size(nbytes)
        if fill is not None:
            if isinstance(fill, bytes):
                self.buffer[self.pos:self.pos+nbytes] = (fill * (nbytes // len(fill) + 1))[0:nbytes]
            else:
                self.buffer[self.pos:self.pos+nbytes] = bytes(nbytes)
        self.pos += nbytes

    def put_data(self, data):
        """
        Just put all the data into record.
        """
        self.assert_bor()
        size = len(data)
        self._check_buf_size(size)
        self.buffer[0:size] = data
        self.pos = size

    def put_n(self, data, dtype = None, order = 'F', offset = None):
        """
        Write numpy object to buffer

        KEYWORDS:
            order - output order of array
                    `None` - use data default
                    default is `F`
            dtype - output data type
        """
        if not isinstance(data, np.ndarray):
            data = np.array(data)
        if dtype is None:
            dtype = data.dtype
        else:
            dtype = np.dtype(dtype)
        if order is None:
            if data.flags.fnc:
                order = 'F'
            else:
                order = 'C'
        assert order in ('F', 'C')
        new_data = np.ndarray(data.shape,
                              dtype = dtype,
                              order = order)
        new_data[()] = data[()]
        data = new_data
        if not data.flags.c_contiguous:
            data = np.ndarray(data.shape,
                              dtype = data.dtype,
                              buffer = data.data,
                              order = 'C')
        if self.swapbyteorder:
            data.byteswap(True)
        nbytes = data.nbytes
        self._check_buf_size(nbytes, offset)
        if offset is not None:
            if offset > 0:
                p = self.pos + offset
            else:
                p = -offset
        else:
            p = self.pos
        self.buffer[p:p+nbytes] = data.data.tobytes()
        if offset is None:
            self.pos += nbytes

    def put_1d(self, data, dtype=np.float64, lead=0, tail=0, offset = None):
        """Write a 1D np array padded with 0 on either side as specified.  Do not write padding.
        """
        if not isinstance(data, np.ndarray):
            data = np.array(data)
        if dtype is None:
            dtype = data.dtype
        else:
            dtype = np.dtype(dtype)
        items = prod(data.shape) - lead - tail
        new_data = np.ndarray(items, dtype = dtype)
        new_data[0:items] = data.flat[lead:lead+items]
        data = new_data
        if self.swapbyteorder:
            data.byteswap(True)
        nbytes = data.nbytes
        self._check_buf_size(nbytes, offset)
        if offset is not None:
            if offset > 0:
                p = self.pos + offset
            else:
                p = -offset
        else:
            p = self.pos
        self.buffer[p:p+nbytes] = data.data.tobytes()
        if offset is None:
            self.pos += nbytes

    @staticmethod
    def get_s_len(s, codec = 'cp437', strip = False):
        """
        Return length of string after encoding.

        If parameter is an array, return array of same shape.
        If parameter is not an np.ndarray, return (nested) list.

        PARAMETERS
        codec - default: 'cp437' used to encode
        """
        t = type(s)
        if not isinstance(s, np.ndarray):
            s = np.array(s, dtype = np.object)
        l = np.ndarray(s.shape, dtype = np.int)
        sflat = s.flat
        lflat = l.flat
        if strip:
            for i in range(len(sflat)):
                lflat[i] = len(sflat[i].strip().encode(codec))
        else:
            for i in range(len(sflat)):
                lflat[i] = len(sflat[i].encode(codec))
        if not issubclass(t, np.ndarray):
            l = l.tolist()
        else:
            l = l[()]
        return l

    def put_s(self, s, length = None, fill = b'\x00', codec = 'cp437', order = 'F', strip = False, offset = None):
        """write string (array) to buffer

        KWARGS
        length - >0: length of string - fill/truncate
                 -1: find and max length
                  None: write actual length of each string
                  np.ndarray: length of strings if match shape
                  (TODO - extend in missing dimesions?)
        fill - pattern (not encoded), memory data if None
        codec - default 'cp437'
        order - of data to write to buffer, default is 'F'
        offset - relative to current location if positive
                 relative to beginning of buffer if negative (abs value)
                 `None` - no offset, advnace buffer
        """
        if order is None:
            if data.flags.fnc:
                order = 'F'
            else:
                order = 'C'
        assert order in ('F', 'C')
        if not isinstance(s, np.ndarray):
            s = np.array(s, dtype = np.object, order = order)
        # create length array
        try:
            if length is None or length == -1:
                l = self.get_s_len(s)
                if length == -1:
                    l = np.max(l)
            else:
                l = length
        except ValueError:
            l = length
        if not isinstance(l, np.ndarray):
            l = np.array(l, dtype = np.int)
        if prod(l.shape) == 1:
            l = np.array(l.flat[0])
        if l.shape == ():
            l = np.tile(l, s.shape)
        if order is 'F' and not s.flags.f_contiguous:
            s = s.copy(order = 'F')
        if order is 'F' and not l.flags.f_contiguous:
            l = l.copy(order = 'F')
        if not s.flags.c_contiguous:
            s = np.ndarray(s.shape,
                           dtype = s.dtype,
                           buffer = s.data,
                           order = 'C')
        if not l.flags.c_contiguous:
            l = np.ndarray(l.shape,
                           dtype = l.dtype,
                           buffer = l.data,
                           order = 'C')
        nbytes = np.sum(l)
        self._check_buf_size(nbytes, offset)
        if offset is not None:
            if offset > 0:
                p = self.pos + offset
            else:
                p = -offset
        else:
            p = self.pos
        if prod(l.shape) > 0:
            lmax = np.max(l)
        else:
            lmax = 0
        f = (fill * (lmax // len(fill) + 1))[:lmax]
        for si, li in zip(s.flat, l.flat):
            d = si.encode(codec)
            n = min(len(d), li)
            self.buffer[p:p+n] = d[:n]
            if fill is not None and n < li:
                self.buffer[p+n:p+li] = f[:li-n]
            p += li
        if offset is None:
            self.pos += nbytes
            assert p == self.pos, "inconsitency in written data"

    # binary data

    def put_buf(self, data, length = None, order = 'F', fill = b'\x00', offset = None):
        """Write array/list of raw data pieces of equal length to buffer.

        ARGS:
        data - array/scalar to be written

        KEYWORDS:
        length - of data pieces, truncate/fill, if `None` use max value
        order  - of junks in written array, default is 'F'
        fill   - default is \\x00

        """
        if fill is not None:
            assert isinstance(fill, bytes) , \
                   "Only bytes-type fill allowed."
        if length is None:
            dtype = np.dtype(np.bytes_)
        else:
            dtype = np.dtype((np.bytes_, length))
        data = np.array(data,
                        dtype = dtype,
                        order = order)
        if not fill in (None, b'\x00'):
            if length is None:
                length = data.dtype.itemsize
            f = (fill * (length // len(fill) + 1)) [:length]
            # array operations for concatenation do not work in numpy 1.11
            d = data.flat
            for i in range(prod(d.shape)):
                d[i] += f
        if not data.flags.c_contiguous:
            data = np.ndarray(data.shape,
                              dtype = data.dtype,
                              buffer = data,
                              order = 'C')
        nbytes = data.nbytes
        self._check_buf_size(nbytes, offset)
        if offset is not None:
            if offset > 0:
                p = self.pos + offset
            else:
                p = -offset
        else:
            p = self.pos
        self.buffer[p:p+nbytes] = data.tobytes()
        if offset is None:
            self.pos += nbytes


    # ========================================
    # application-specific routines
    # ========================================

    def put_kep_parm(self, data):
        """Write a kepler parameter binary list with 32 bit integers."""
        count = len(data)
        value = np.zeros(
            count,
            dtype=np.float64)
        ivalue = np.ndarray(
            count,
            buffer=value.data.cast('b'),
            offset=4,
            dtype=np.int32,
            strides=8)
        for i,d in enumerate(data):
            if d.dtype == np.int32:
                ivalue[i] = d
            else:
                value[i] = d
        if self.swapbyteorder:
            value.byteswap(True)
        p = self.pos
        nbytes = value.nbytes
        self._check_buf_size(nbytes)
        self.buffer[p:p+nbytes] = value.data.tobytes()
        self.pos += nbytes

    def put_kep_parm64(self, data):
        """Write a kepler parameter binary list with 64 bit integers."""
        count = len(data)
        if count == 0:
            return
        value = np.zeros(
            count,
            dtype=np.float64)
        ivalue = np.ndarray(
            count,
            buffer=value.data.cast('b'),
            dtype=np.int64)
        for i,d in enumerate(data):
            if d.dtype == np.int64:
                ivalue[i] = d
            else:
                value[i] = d
        if self.swapbyteorder:
            value.byteswap(True)
        p = self.pos
        nbytes = value.nbytes
        self._check_buf_size(nbytes)
        self.buffer[p:p+nbytes] = value.data.tobytes()
        self.pos += nbytes

    def put_f8_kep_i4(self, data):
        """Write i4 in f8 array for kepler.

        Pass the f8 dimension.

        Half the space seems wasted the way KEPLER treats this, the
        entire second half of each array is empty.

        Here we shall just fill up the 2nd part of the array and
        write the passed dimension.

        Byteswap is only needed on i4 level (see read routine).
        """
        self.put_n(data)
        self.skip_bytes(data.nbytes, fill=b'\x00')

    # dummy routine to allow data IO code below

    def write(self):
        """
        Provide interface for writing data to whereever.
        """
        raise NotImplementedError("Writing Data not implemented.")

    # =======================================================================
    # data IO
    # =======================================================================

    def write_bytes(self, *args, **kwargs):
        """Write numpy empty bytes to file"""
        self.assert_bor()
        self.skip_bytes(*args, **kwargs)
        self.write()

    def write_data(self, *args, **kwargs):
        """Write plain buffer to file"""
        self.assert_bor()
        self.put_data(*args, **kwargs)
        self.write()

    def write_n(self, *args, **kwargs):
        """Write numpy scalar/array to file"""
        self.assert_bor()
        kwargs['offset'] = None
        self.put_n(*args, **kwargs)
        self.write()

    def write_1d(self, *args, **kwargs):
        """Write 1d padded numpy scalar/array to file"""
        self.assert_bor()
        kwargs['offset'] = None
        self.put_1d(*args, **kwargs)
        self.write()

    def write_s(self, *args, **kwargs):
        """Write (numpy) string (array) to file"""
        self.assert_bor()
        kwargs['offset'] = None
        self.put_s(*args, **kwargs)
        self.write()

    def write_buf(self, *args, **kwargs):
        """Write array/list of raw data pieces to file"""
        self.assert_bor()
        kwargs['offset'] = None
        self.put_buf(*args, **kwargs)
        self.write()

    # application-specific routines

    def write_kep_parm(self, data):
        """
        write kepler parm array to file
        """
        self.assert_bor()
        self.put_kep_parm(data)
        self.write()

    def write_f8_kep_i4(self, data):
        """Write i4 in f8 array for kepler.

        Pass the f8 dimension.

        Half the space seems wasted the way KEPLER treats this, the
        entire second half of each arry is empty.

        Here we shall just fill up the 2nd part of the array and
        write the passed dimension.

        Byteswap is only needed on i4 level (see read routine).
        """
        self.assert_bor()
        self.put_f8_kep_i4(data)
        self.write()

def _f_store_n(p = None, dt = None, **_kwargs):
    def _f(self, *args, **kwargs):
        kwargs['dtype'] = dt
        kwargs.setdefault('order', 'F')
        p(self, *args, **kwargs)
    _f.dt = dt
    _f.p = p
    return _f

def _f_store_n1d(p = None, dt = None, **_kwargs):
    def _f(self, data, *args, **kwargs):
        kwargs['dtype'] = dt
        p(self, data, *args, **kwargs)
    _f.p = p
    _f.dt = dt
    return _f

def _f_store_n1d_(p = None, dt = None, lead = 0, tail = 0, **_kwargs):
    def _f(self, data, *args, **kwargs):
        kwargs['dtype'] = dt
        kwargs['lead'] = lead
        kwargs['tail'] = tail
        p(self, data, *args, **kwargs)
    _f.p = p
    _f.dt = dt
    _f.lead = lead
    _f.tail = tail
    return _f


for t in _np_types:
    kw = dict(
        cls = DataOutputBuffer,
        t = t)
    _set_method(
        fn = _f_store_n,
        parent = 'put_n',
        name = 'put_{t}',
        doc = """Write numpy {dn} to buffer at offset relative to
        current position.\n
        Does not advance buffer pointer.""",
        **kw)
    _set_method(
        fn = _f_store_n,
        parent = 'write_n',
        name = 'write_{t}',
        doc = "Write numpy {dn} array to file as record.",
        **kw)

    _set_method(
        fn = _f_store_n1d,
        parent = 'put_1d',
        name = 'put_{t}_1d',
        doc = """Write a 1D numpy {dn} array padded with 0 as specified to buffer.  Padding is not written.""",
        **kw)
    _set_method(
        fn = _f_store_n1d,
        parent = 'write_1d',
        name = 'write_{t}_1d',
        doc = """Write a 1D numpy {dn} array padded with 0 as specified to file as record.  Padding is not written.""",
        **kw)

    _set_method(
        fn = _f_store_n1d_,
        parent = 'put_1d',
        name = 'put_{t}_1d_0',
        doc = """Write a 1D numpy {dn} array padded with one element at beginning.  Padding is not written.""",
        extra_kw = dict(lead=1, tail=0),
        **kw)
    _set_method(
        fn = _f_store_n1d_,
        parent = 'write_1d',
        name = 'write_{t}_1d_0',
        doc = """Write a 1D numpy {dn} array padded with one element at beginning to file as record.  Padding is not written.""",
        extra_kw = dict(lead=1, tail=0),
        **kw)

    _set_method(
        fn = _f_store_n1d_,
        parent = 'put_1d',
        name = 'put_{t}_1d_n',
        doc = """Write a 1D numpy {dn} array padded with one element at end.  Padding is not written.""",
        extra_kw = dict(lead=0, tail=1),
        **kw)
    _set_method(
        fn = _f_store_n1d_,
        parent = 'write_1d',
        name = 'write_{t}_1d_n',
        doc = """Write a 1D numpy {dn} array padded with one element at end to file as record.  Padding is not written.""",
        extra_kw = dict(lead=0, tail=1),
        **kw)

    _set_method(
        fn = _f_store_n1d_,
        parent = 'put_1d',
        name = 'put_{t}_1d_0n',
        doc = """Write a 1D numpy {dn} array padded with one element at begiining and end each.  Padding is not written.""",
        extra_kw = dict(lead=1, tail=1),
        **kw)
    _set_method(
        fn = _f_store_n1d_,
        parent = 'write_1d',
        name = 'write_{t}_1d_0n',
        doc = """Write a 1D numpy {dn} array padded with one element at beginning and end each to file as record.  Padding is not written.""",
        extra_kw = dict(lead=1, tail=1),
        **kw)

#=======================================================================

class DataWriter(DataOutputBuffer):
    """
    Class for writing 'unformatted' binary files.

    File names ending with .gz, .xz, .bz2 will be automatically
    compressed.

    For .gz this may fail, however, if the file is bigger than 2GB or
    4GB.
    """
    # TODO: add (file) truncate functionallity

    def __init__(self,
                 filename, *args, **kwargs):
        """
        Initialize data fields and open file.

        Optionally the byte order can be specified.
        The default is big endian.
        """

        # TODO: Add append mode.
        # not sure how/whether this will work with compressed files

        super().__init__(*args, **kwargs)
        self.open(filename)

    def open(self, filename):
        """
        Open the file for writing.
        """
        self.filename = os.path.expandvars(os.path.expanduser(filename))

        if self.filename.endswith('.gz'):
            self.compressed = True
            self.compress_mode = 'gz'
            self.file = gzip.open(filename,'wb')
        elif self.filename.endswith('.bz2'):
            self.compressed = True
            self.compress_mode = 'bz2'
            self.file = bz2.BZ2File(filename,'wb',2**16)
        elif self.filename.endswith('.xz'):
            self.compressed = True
            self.compress_mode = 'xz'
            self.file = LZMAFile(self.filename,'wb')
        else:
            self.file = open(filename,'wb',-1)
            self.compressed = False
            self.compress_mode = None
        self._init()

    def _init(self):
        """Initialize the file position and data to empty."""
        super()._init()
        self.fpos = 0

    def close(self):
        """Close the file."""
        if self.pos != 0:
            self.write()
        self.file.close()

    def rewind(self):
        """Rewind the file."""
        self.file.seek(0, os.SEEK_SET)
        self._init()

    def write(self):
        """
        Write the data record to file.
        """
        self._write()

    def _write(self):
        """
        Write a data record to file.
        """
        self.file.write(self.buffer[0:self.pos])

        self.fpos += self.pos
        self.pos = 0

    # context manager interface

    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is None:
            self.close()
        return False

#=======================================================================

class FortranWriter(DataWriter, FortranSpecs):
    """
    Class for writing 'unformatted' Fortran binary files.

    Based on DataWriter, automatic compression support.
    """

    # I still need to find out how to efficiently work with a buffer
    # for writing internally - how to extend it, etc.  The plan is to
    # always write an entire record at once.  Currently it would
    # appear the largest FORTRAN files I have come with records for
    # ppnb in KEPLER, 5000 isotopes * 2000 zones * 8 bytes ... 80 MB
    # if pushing it ... usually < 16 MB.  So we could start with that,
    # then extend if ever needed.

    # TODO - add functionallity for record-based backward skipping
    # (potentially need to tuncate file on close)

    def __init__(self, *args, reclen = 4, **kwargs):
        self._set_reclen(reclen)
        super().__init__(*args, **kwargs)

    def _init(self):
        """Initialize the file position and data to empty."""
        super()._init()
        self.rpos = 0

    def _set_byteorder(self, *args, **kwargs):
        super()._set_byteorder(*args, **kwargs)
        self.reclen_dtype = self.reclen_dtype.newbyteorder(self.byteorder)

    def _write(self):
        """
        Write a data record to file.
        """
        self._write_reclen()
        self.file.write(self.buffer[0:self.pos])
        self._write_reclen()

        self.fpos += self.pos + 2 * self.fortran_reclen
        self.rpos += 1
        self.pos = 0

    def _write_reclen(self):
        """Write the record length."""
        self.file.write(np.array(
            self.pos,
            dtype = self.reclen_dtype).data.tobytes())
