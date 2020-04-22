"""
Classes for reading UNIX unformatted FORTRAN files.
"""

# TODO
# * for general reading, load needs to be able to specify record size?

import os
import sys
import gzip
import bz2
import lzma
import glob
import itertools

import numpy as np

from types import MethodType
from copy import copy

from .utils import prod, cumsum0_enum_range_iter, xz_file_size

from .errors import RecordSizeError, ReadError, StringReadError, FileReadError, \
     DataFileError, RecLenError

from .common import _np_types, _type_map, _set_method

from .defaults import FortranSpecs

class DataInputBuffer(object):
    """
    Provide basic IO for reading from data buffer.

    In the routine descriptions 'record length' referse to the size oc
    the current data buffer.

    The routines allow to specify byte order in the buffer.  This is
    used when readling data from files that were written with
    different byte order.

    TODO - specify naming conventions for routines.
    """

    default_byte_order = '='

    sys_is_le = sys.byteorder == 'little'
    native_byteorder = '<' if sys_is_le else '>'

    def __init__(self,
                 byteorder = None,
                 **kwargs):
        self.reclen = 0
        self.pos = 0
        self.data = b''
        self._set_byteorder(byteorder = byteorder)

    def _set_byteorder(self,
                       byteorder = None):
        """
        set up all data types for deserved byte order
        """
        if byteorder is None:
            byteorder = self.default_byte_order
        if byteorder == "=":
            byteorder = self.native_byteorder
        self.swapbyteorder = byteorder != self.native_byteorder
        self.byteorder = byteorder

    def _reverse_byteorder(self):
        """
        Reverse data type byte order
        """
        byteorder = '<' if self.byteorder == '>' else '>'
        self._set_byteorder(byteorder)

    def eor(self):
        """Return whether position is end of current record."""
        return self.pos == self.reclen

    def assert_eor(self):
        """Throw exception if current position is not end of record.

        This can be used to deterim whether all data has been read,
        i.e., as a consistency check of the read data sizes"""
        if not self.eor():
            raise RecordSizeError(self.reclen, self.pos)

    def check_buf(self, length, offset=0):
        """Check whether length bytes are still available on buffer.

        This can can also be overwritten to load more data in case the
        record is not loaded at once.

        Negative offsets are relaitive to end of buffer.
        """
        offset = int(offset)
        if length < 0:
            raise ReadError(
                '[check_buf] Reading length must not be negative: length = {}.'.format(
                 length))
        if offset < 0:
            pos = self.reclen + offset
        else:
            pos = self.pos + offset
        if pos < 0:
            raise ReadError(
                '[check_buf] Reading before beginning of record by {} bytes.'.format(
                 - (offset + self.pos)))
        if pos + length > self.reclen:
            raise ReadError(
                '[check_buf] Reading past end of record by {} bytes.'.format(
                offset + length + self.pos - self.reclen))

    def get_buf_avail(self, size=1, offset=0, tail=0, return_rest = True):
        """
        detemine number of available units
        """
        if size <= 0:
            raise ReadError('[buf_avail] size must be > 0')
        if offset < 0:
            pos = self.reclen + offset
        else:
            pos = self.pos + offset
        if tail is None:
            tail = 0
        if tail < 0:
            raise ReadError('[buf_avail] tail must be >= 0')
        nbytes = self.reclen - pos - tail
        if return_rest:
            return divmod(nbytes, size)
        return nbytes // size

    # dummy interface routine for reading
    def load(self, *args, **kwargs):
        """Fill buffer with new data"""
        raise NotImplementedError()

    def skip_bytes(self, length=None):
        """Skip number of bytes on read."""
        if length is None:
            length = self.reclen - self.pos
        self.check_buf(length)
        self.pos += length

    def resolve_free_dimenions(self, dim, dtype, offset, tail):
        """deal with free dimensions

        check space.

        return
           (dim, length)

           dim - resolved dimensions
           length - length of structure in bytes

        TODO - treatment of None, np.newaxis?
        """
        if dim is None:
            dim = ()
        try:
            dim = tuple(dim)
        except TypeError:
            if np.isscalar(dim):
                dim = (int(dim),)
            else:
                raise ReadError('[dimensions] Cannot interpret dim argument.')

        length = prod(dim) * dtype.itemsize
        nfree = dim.count(-1)
        if nfree > 1:
            raise ReadError('[dimensions] More than one free dimension.')
        elif nfree == 1:
            ifree = dim.index(-1)
            recsize = -length
            if tail is None:
                ntail = 0
            else:
                ntail = tail
            nrec, bytes_remain = get_buf_avail(
                size=recsize,
                offset=offset,
                tail=ntail)
            if tail is not None and bytes_remain != tail:
                raise ReadError('[dimensions] Cannot match free dimension.')
            length = nrec * recsize
            dim = dim[:ifree] + (nrec,) + dim[ifree+1:]
        else:
            self.check_buf(length, offset)
        return dim, length

    ############################################################
    # peak routines - fundamental loading
    ############################################################

    def peek_data(self):
        """Just get all of the record data."""
        if self.pos > 0:
            raise self.ReadError('[get_data] Not at beginning of record')
        return self.peek_buf()

    def peek_buf(self, dim=(), length=None, offset=0, order='F', output='C', truncate=False, tail=None):
        """Read some raw data.

        Default length is rest of buffer unless non-scalar dimesnion
        is specified, in which case length needs to be specified.

        KEYWORDS
        output - `list` return list rather then numpy array
        order  - order of data in buffer/file
        offset - offset in buffer/file
        truncate - allow numpy to truncate 0 bytes at end of buffer
                   (seems to be not what we usually want)
                   will create bytes object for each array entry
        """
        if length is None:
            assert dim == (), '[{}] dimension needs to be specified'.format(self.__class__.__name__)
            length = self.reclen - self.pos - offset

        dtype = np.dtype((np.bytes_, length))

        dim, _ = self.resolve_free_dimenions(dim, dtype, offset, tail)

        value = np.ndarray(
            dim,
            dtype = dtype,
            offset = self.pos + offset,
            buffer = self.data,
            order = order)

        # np.bytes truncates buffer items on access, we need to create objects
        if not truncate:
            value_bytes = value
            value = np.ndarray(dim, dtype = np.object)
            bflat = value_bytes.flat
            vflat = value.flat
            for i, b in enumerate(bflat):
                vflat[i] = b + b'\0' * (length - len(b))

        if output == 'list':
            return value.copy().tolist()
        if output is None:
            output = order
        value = value.copy(order = output)
        return value[()]

    def peek_n(self, dim=None, dtype=None, *, offset=0, order='F', output='C', tail=None):
        """Get one numpy value at location offset relative to current position.  Do not advance buffer pointer.

        Arguments
          dim and type may be switched.

        Keyword Arguments
          dim - tuple of dimensions, may contain one '-1' value, default is ()
          dtype - numpy.dtype to use, default is  np.float64

        KW only
          offset - offset in bytes relative to current location
          order  - order of data on disk
          output - order of data array to be returned
                   `None` defaults to `order`
                   `'list'`
                   default is `'C'`
          tail - number of bytes not to read at end of record;
                 used when one of the dimensions is -1
                 None: ignore
        """
        if type(dim) == type:
            dtype, dim = dim, dtype
        if dtype is None:
            dtype = np.float64
        if not isinstance(dtype, np.dtype):
            dtype = np.dtype(dtype)

        dim, _ = self.resolve_free_dimenions(dim, dtype, offset, tail)

        value = np.ndarray(
            dim,
            buffer = self.data,
            offset = self.pos + offset,
            dtype = dtype,
            order = order)
        if output is None:
            output = order
        value = value.copy(order=output)
        if self.swapbyteorder and dtype.byteorder == '=':
            value.byteswap(True)
        if output == 'list':
            return value.tolist()
        return value[()]

    # fundamental routine for 1D peaks with padding
    # TODO - generalise to multi-D padding
    #        pad = ((dim0pad0, dim0padn), (dim1pad0, dim1padn), ...)
    def peek_1d(self, length, dtype=np.float64, lead=0, tail=0, output = None):
        """Read a 1D np array and pad with 0 on either side as specified.  Do not advance buffer pointer."""
        if not isinstance(dtype, np.dtype):
            dtype = np.dtype(dtype)
        itemsize = dtype.itemsize
        nb = itemsize * length
        self.check_buf(nb)
        totlen = length + lead + tail
        value = np.ndarray(totlen, dtype=dtype)
        p = self.pos
        i0 = itemsize * lead
        value.view(np.uint8).data[i0:i0+nb] = self.data[p:p+nb]
        if self.swapbyteorder:
            value.byteswap(True)
        value[0:lead]=0
        value[totlen-tail:totlen]=0
        if output == 'list':
            return value.tolist()
        return value

    # STRING FUNCTIONS

    def peek_sln(self, *args, offset = 0, codec = 'cp437', order='F', output='C', **kwargs):
        """
        Read an array of stings of varying lengths.

        Parameter 'lengths' is required, as first positional argument
        or kw argument.  'lengths' are the byte sizes on disk, not the
        decoded unicode byte lengths.

        Return a ndarray of type object filled with strings.

        The maximum unit ode string size may depend on codec and
        content and cannot be determined a priory from the byte size on file.

        ARGS
        [[dim], length]
        dim - dimension, default is ()
        length - of string, default is 1

        KEYWORDS
        order  - order of input array
        output - order of output array; `'list'` to return (nested) list
        codec  - `'cp437'`
        offset - relative stars in buffer, 0
        """
        if 'lengths' in kwargs:
            lengths = kwargs.pop('lengths')
        elif len(args) == 0:
            raise StringReadError("Parameter 'lengths' is required.")
        else:
            lengths = args[0]
            args = args[1:]
        if len(args) > 0:
            raise TypeError('Wrong number of arguments supplied.')
        if len(kwargs) > 0:
            raise TypeError('Invalid keyword arguments supplied: {}.'.format(kwargs))

        if not isinstance(lengths, np.ndarray):
            lengths = np.array(lengths, dtype = np.int)
        value = np.ndarray(lengths.shape, dtype = np.object, order = order)
        flat = value.flat
        nbytes = np.sum(lengths)
        self.check_buf(nbytes, offset)
        for k, i, j in cumsum0_enum_range_iter(lengths.flat, self.pos + offset):
            flat[k] = self.data[i:j].decode(codec)
        if not value.flags.c_contiguous:
            value = np.ndarray(value.shape,
                               dtype = value.dtype,
                               buffer = value,
                               order = 'C')
        if output == 'list':
            return value.copy().tolist()
        if output is None:
            output = order
        value = value.copy(order = output)
        return value[()]

    def peek_s(self, *args, **kwargs):
        """Same as peek_sn, but return (nested) list."""
        return self.peek_sn(*args, output = 'list', **kwargs)

    def peek_sn(self, *args, strip = False, offset = 0, order = 'F', output = 'C', codec = 'cp437', tail = None, **kwargs):
        """Read srings of fixed length (in buffer).

        call signature is peak_s([[dim,] length], ...)

        'dim' and 'length' may also be provided as keyword arguments.

        RETURNS
          numpy ndarray of unicode elelemts of requested length

        ARGS/KWARGS
        dim - dimension of array, default: () for scalars
        length - of string in bytes to read from file, default: 1
        offset - where to read, default: 0
        codec - codec for conversion from bytes to unicode, default cp437
        strip - whether to strip string, default: False
        output - 'C', 'F', or 'list' - memory order for multi-D
        """
        if 'length' in kwargs:
            length = kwargs.pop('length')
        else:
            try:
                length = args[-1]
                args = args[:-1]
            except:
                length = 1
        if 'dim' in kwargs:
            dim = kwargs.pop('dim')
        else:
            if len(args) > 0:
                dim = args[0]
                args = args[1:]
            else:
                dim = ()
        if dim is None:
            dim = ()
        try:
            dim = tuple(dim)
        except TypeError:
            dim = (dim,)
        if len(args) > 0:
            raise TypeError('Wrong number of arguments supplied.')
        if len(kwargs) > 0:
            raise TypeError('Invalid keword arguments supplied: {}.'.format(kwargs))

        dtype = np.dtype((np.bytes_, length))

        dim, _ = self.resolve_free_dimenions(dim, dtype, offset, tail)

        value_bytes = np.ndarray(
            dim,
            buffer = self.data,
            offset = self.pos,
            dtype = dtype,
            order = order)
        value = np.ndarray(dim, dtype = np.dtype((np.str, length)))
        bflat = value_bytes.flat
        vflat = value.flat
        for i, v in enumerate(bflat):
            vflat[i] = v.decode(codec)
            if strip is True:
                vflat[i] = vflat[i].strip()
        if output == 'list':
            return value.copy().tolist()
        if output is None:
            output = order
        value = value.copy(order = output)
        return value[()]

    def peek_s(self, *args, **kwargs):
        """Same as peek_sn, but return (nested) list."""
        return self.peek_sn(*args, output = 'list', **kwargs)

    ############################################################
    # get routines - loading and advancing buffer
    ############################################################

    def get_data(self):
        """Just get all of the record data."""
        data = self.peek_data()
        self.pos += len(data)
        return data

    def get_buf(self, dim=(), length=None, **kwargs):
        """Read some raw data."""
        value = self.peek_buf(dim=dim, length=length, offset=0, **kwargs)
        if dim == ():
            self.pos += len(value)
        else:
            self.pos += prod(dim) * length
        return value

    def get_data(self):
        """Just get all of the record data."""
        value = self.peek_data()
        self.pos = self.reclen
        return value

    # this one should be able to replace all the multiple-dim routines.
    def get_n(self, *args, **kwargs):
        """Read an np array of type dtype.
        """
        value = self.peek_n(*args, **kwargs)
        self.pos += value.nbytes
        return value
    get_n.__doc__ += '\n'.join(peek_n.__doc__.splitlines()[1:])

    # fundamental routine for 1D loads with padding
    def get_1d(self, length, dtype=np.float64, **kwargs):
        """Read a 1D np array and pad with 0 on either side as specified."""
        if not isinstance(dtype, np.dtype):
            dtype = np.dtype(dtype)
        itemsize = dtype.itemsize
        nb = itemsize * length
        value = self.peek_1d(length, dtype=dtype, **kwargs)
        self.pos += nb
        return value

     # string get routines

    def get_sln(self, *args, **kwargs):
        """
        Read an array of stings of varying lengths.

        Here we return a ndarray of type object filled with strings.
        """
        lengths = kwargs.get('lengths', args[0])
        kwargs['offset'] = 0
        value = self.peek_sln(*args, **kwargs)
        self.pos += np.sum(lengths)
        return value

    def get_sl(self, *args, **kwargs):
        """Same as get_sln, but return (nested) list."""
        return self.get_sln(*args, output='list', **kwargs)

    def get_sn(self, *args, **kwargs):
        """Similar to peek_sn, but read data from current position, no offset allowed."""
        # ensure we do not pass an offset
        kwargs['offset'] = 0
        value = self.peek_sn(*args, **kwargs)
        try:
            length = kwargs.get('lenght', args[-1])
        except:
            length = 1
        nbytes = value.size * length
        self.pos += nbytes
        return value

    def get_s(self, *args, **kwargs):
        """Same as get_sn, but return (nested) list."""
        return self.get_sn(*args, **kwargs).tolist()

    ############################################################
    # direct IO routines
    ############################################################

    def load_1d(self, *args, **kwargs):
        """Load a 1D np array and pad with 0 on either side as specified."""
        self.load()
        value = self.get_1d(*args, **kwargs)
        self.assert_eor()
        return value

    def load_data(self):
        """Just load all of the record data."""
        self.load()
        data = self.get_data()
        self.assert_eor()
        return data

    def load_buf(self, *args, **kwargs):
        """Loads some raw data."""
        self.load()
        value = self.get_buf(*args, **kwargs)
        self.assert_eor()
        return value

    def load_n(self, *args, **kwargs):
        """Loads some np array."""
        self.load()
        value = self.get_n(*args, **kwargs)
        self.assert_eor()
        return value

    # string routines

    def load_sln(self, *args, **kwargs):
        """Load and read string array of variable length and return as numpy array."""
        self.load()
        value = self.get_sln(*args, **kwargs)
        self.assert_eor()
        return value

    def load_sl(self, *args, **kwargs):
        """Load and read string array of variable length and return as numpy array."""
        return self.load_sln(*args, **kwargs).tolist()

    def load_sn(self, *args, **kwargs):
        """Load and read string array of fixed length and return as numpy array."""
        self.load()
        value = self.get_sn(*args, **kwargs)
        self.assert_eor()
        return value

    def load_s(self, *args, **kwargs):
        """Load and read string array of fixed length and return as numpy array."""
        return self.load_sn(*args, **kwargs).tolist()

    ############################################################
    # the following are application-specific routines but could serve
    # as an example for other applications.  If they were moved into
    # derived class, iterrecords would no longer work.
    ############################################################

    def get_kep_parm(self,list):
        """Read a kepler parameter binary list."""
        count = len(list)
        value = np.ndarray(
            count,
            buffer=self.data,
            offset=self.pos,
            dtype=np.float64)
        value = value.copy()
        if self.swapbyteorder:
            value.byteswap(True)
        ivalue = np.ndarray(
            count,
            buffer=value,
            offset=4,
            dtype=np.int32,
            strides=8)
        self.pos += 8 * count
        return (value[i] if l == 1 else ivalue[i] for i,l in enumerate(list))

    def get_kep_parm64(self, list):
        """Read a 64-bit kepler parameter binary list."""
        count = len(list)
        value = np.ndarray(
            count,
            buffer=self.data,
            offset=self.pos,
            dtype=np.float64)
        value = value.copy()
        if self.swapbyteorder:
            value.byteswap(True)
        ivalue = np.ndarray(
            count,
            buffer=value.data,
            dtype=np.int64)
        self.pos += 8 * count
        return (value[i] if l == 1 else ivalue[i] for i,l in enumerate(list))

    def get_f8_kep_i4n(self,dim):
        """Read i4 in f8 array for kepler.

        Pass the f8 dimension.

        Half the space seems wasted the way KEPLER treats this, the
        entire second half of each arry is empty.

        Here we shall just discard the 2nd part of the array and only
        return the requested dimension.

        Usually one would first read the f8, do a byteswap, then use
        the buffer for the integers, however, KEPLER stores the i4 in
        the f8 in the same format a big endian system would have, and
        hence the byte-swap is only needed on the i4 level.
        """
        value = np.ndarray(
            dim,
            buffer = self.data,
            offset = self.pos,
            dtype = np.int32,
            order = 'F')
        value = value.copy()
        if self.swapbyteorder:
            value.byteswap(True)
        self.pos += 8 * prod(dim)
        return value

    # application-specfic load routines

    def load_f8_kep_i4n(self, dim):
        """Load and read i4 in f8 array for kepler."""
        self.load()
        value = self.get_f8_kep_i4n(dim)
        self.assert_eor()
        return value

# generate convenience functions

# p stands for 'parent' function
# dt  for 'data type'

def _f_peek_n(p = None, dt = None, **_kwargs):
    def _f(self, *args, **kwargs):
        kwargs['dtype'] = dt
        kwargs.setdefault('order', 'F')
        return p(self, *args, **kwargs)
    _f.p = p
    _f.dt = dt
    return _f

def _f_peek(p = None, **kwargs):
    def _f(self, dim = (), offset = 0):
        return p(self, dim=dim, offset=offset).tolist()
    _f.p = p
    return _f

def _f_get_n(p = None, **_kwargs):
    def _f(self, dim = ()):
        value = p(self, dim=dim)
        self.pos += value.nbytes
        return value
    _f.p = p
    return _f

def _f_get(p = None, **_kwargs):
    def _f(self, dim = ()):
        return p(self, dim=dim).tolist()
    _f.p = p
    return _f

def _f_load_n(p = None, **_kwargs):
    def _f(self, dim = ()):
        self.load()
        value = p(self, dim=dim)
        self.assert_eor()
        return value
    _f.p = p
    return _f

def _f_load(p = None, **_kwargs):
    def _f(self, dim = ()):
        return p(self, dim=dim).tolist()
    _f.p = p
    return _f

def _f_check_illeagal_argument(argnames, kwargs):
    if not isinstance(argnames, str):
        argnames = (argnames,)
    for argname in argnames:
        assert argname not in kwargs, 'illegal argument {argname}: "{argval}"' \
               .format(argname=argname, argval=kwargs[argname])

def _f_peak_n1d(dt = None, **_kwargs):
    def _f(self, length, **kwargs):
        _f_check_illeagal_argument('dtype', kwargs)
        return self.peek_1d(length, dtype=dt, **kwargs)
    _f.dt = dt
    return _f

def _f_get_n1d(dt = None, **_kwargs):
    def _f(self, length, **kwargs):
        _f_check_illeagal_argument('dtype', kwargs)
        return self.get_1d(length, dtype=dt, **kwargs)
    _f.dt = dt
    return _f

def _f_load_n1d(dt = None, **_kwargs):
    def _f(self, length, **kwargs):
        _f_check_illeagal_argument('dtype', kwargs)
        return self.load_1d(length, dtype=dt, **kwargs)
    _f.dt = dt
    return _f

def _f_peak_n1d_(dt = None, lead=0, tail=0, **_kwargs):
    def _f(self, length, **kwargs):
        _f_check_illeagal_argument(('dtype', 'lead', 'tail'), kwargs)
        return self.peek_1d(length, dtype=dt, lead=lead, tail=tail, **kwargs)
    _f.dt = dt
    _f.lead = lead
    _f.tail = tail
    return _f

def _f_get_n1d_(dt = None, lead=0, tail=0, **_kwargs):
    def _f(self, length, **kwargs):
        _f_check_illeagal_argument(('dtype', 'lead', 'tail'), kwargs)
        return self.get_1d(length, dtype=dt, lead=lead, tail=tail, **kwargs)
    _f.dt = dt
    _f.lead = lead
    _f.tail = tail
    return _f

def _f_load_n1d_(dt = None, lead=0, tail=0, **_kwargs):
    def _f(self, length, **kwargs):
        _f_check_illeagal_argument(('dtype', 'lead', 'tail'), kwargs)
        return self.load_1d(length, dtype=dt, lead=lead, tail=tail, **kwargs)
    _f.dt = dt
    _f.lead = lead
    _f.tail = tail
    return _f

for t in _np_types:
    kw = dict(
        cls = DataInputBuffer,
        t = t)
    _set_method(
        fn = _f_peek_n,
        parent = 'peek_n',
        name = 'peek_{t}n',
        doc = """Get numpy {dn} array at location offset relative to
        current position.\n
        Does not advance buffer pointer.
        dimesion '()' returns scalar.""",
        **kw)
    _set_method(
        fn = _f_get_n,
        name = 'get_{t}n',
        parent = 'peek_{t}n',
        doc = "Read numpy {dn} array.",
        **kw)
    _set_method(
        fn = _f_get,
        name = 'get_{t}',
        parent = 'get_{t}n',
        doc = "Get one Python {pt} from {dn}.",
        **kw)
    _set_method(
        fn = _f_peek,
        name = 'peek_{t}',
        parent = 'peek_{t}n',
        doc = """Read one numpy {t} at location offset relative to
        current position.  Do not advance buffer pointer.  Return
        Python {pt}""",
        **kw)
    _set_method(
        fn = _f_load_n,
        name = 'load_{t}n',
        parent = 'get_{t}n',
        doc = "Load and read numpy {dn} array.",
        **kw)
    _set_method(
        fn = _f_load,
        name = 'load_{t}',
        parent = 'load_{t}n',
        doc = "Load and read numpy {dn} and return Python {pt}.",
        **kw)

    _set_method(
        fn = _f_peak_n1d,
        name = 'peek_{t}n1d',
        doc = """Read a 1D numpy {dn} array and pad with 0 as specified.  Do not advance buffer pointer.""",
        **kw)
    _set_method(
        fn = _f_get_n1d,
        name = 'get_{t}n1d',
        doc = """Read a 1D numpy {dn} array and pad with 0 as specified.""",
        extra_kw = dict(lead=0, tail=0),
        **kw)
    _set_method(
        fn = _f_load_n1d,
        name = 'load_{t}n1d',
        doc = """Load a 1D numpy {dn} array and pad with 0 as specified.""",
        extra_kw = dict(lead=0, tail=0),
        **kw)

    _set_method(
        fn = _f_peak_n1d_,
        name = 'peek_{t}n1d0',
        doc = """Read a 1D numpy {dn} array and front-pad with one 0.  Do not advance buffer pointer.""",
        extra_kw = dict(lead=1, tail=0),
        **kw)
    _set_method(
        fn = _f_peak_n1d_,
        name = 'peek_{t}n1dn',
        doc = """Read a 1D numpy {dn} array and pad with one 0 at end.  Do not advance buffer pointer.""",
        extra_kw = dict(lead=0, tail=1),
        **kw)
    _set_method(
        fn = _f_peak_n1d_,
        name = 'peek_{t}n1d0n',
        doc = """Read a 1D numpy {dn} array and pad with one 0 on both sides.  Do not advance buffer pointer.""",
        extra_kw = dict(lead=1, tail=1),
        **kw)

    _set_method(
        fn = _f_get_n1d_,
        name = 'get_{t}n1d0',
        doc = """Read a 1D numpy {dn} array and front-pad with one 0.""",
        extra_kw = dict(lead=1, tail=0),
        **kw)
    _set_method(
        fn = _f_get_n1d_,
        name = 'get_{t}n1dn',
        doc = """Read a 1D numpy {dn} array and pad with one 0 at end.""",
        extra_kw = dict(lead=0, tail=1),
        **kw)
    _set_method(
        fn = _f_get_n1d_,
        name = 'get_{t}n1d0n',
        doc = """Read a 1D numpy {dn} array and pad with one 0 on both sides.""",
        extra_kw = dict(lead=1, tail=1),
        **kw)

    _set_method(
        fn = _f_load_n1d_,
        name = 'load_{t}n1d0',
        doc = """Load a 1D numpy {dn} array and front-pad with one 0.""",
        extra_kw = dict(lead=1, tail=0),
        **kw)
    _set_method(
        fn = _f_load_n1d_,
        name = 'load_{t}n1dn',
        doc = """Load a 1D numpy {dn} array and pad with one 0 at end.""",
        extra_kw = dict(lead=0, tail=1),
        **kw)
    _set_method(
        fn = _f_load_n1d_,
        name = 'load_{t}n1d0n',
        doc = """Load a 1D numpy {dn} array and pad with one 0 on both sides.""",
        extra_kw = dict(lead=1, tail=1),
        **kw)

#=======================================================================

class BufferReader(DataInputBuffer):
    """
    Provide basic IO for reading from buffered data.
    """

    def __init__(self,
                 data,
                 byteorder = None,
                 *args,
                 **kwargs):
        """
        Initialize data fields and open file.

        byteorder:
            '<', '>', '=':
                little, big, native endian
                x86 have native '<'
                risc have native '>'
            None:
                try native, then check if possible
                (default)
        """
        self._set_byteorder(byteorder = byteorder)
        self.data = data
        self.reclen = len(data)
        self.pos = 0

    def rewind(self):
        """Rewind the file."""
        pos = 0

    # context manager interface

    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_val, exc_tb):
        return False

#=======================================================================

class DataReader(DataInputBuffer):
    """
    Provide basic IO for reading from buffered data file.
    """

    def __init__(self,
                 filename = None,
                 extension = None,
                 tolerant = False,
                 byteorder = None,
                 *args,
                 **kwargs):
        """
        Initialize data fields and open file.

        byteorder:
            '<', '>', '=':
                little, big, native endian
                x86 have native '<'
                risc have native '>'
            None:
                try native, then check if possible
                (default)

        tolerant:
            True:
                allow partial read of record
            False:
                throw exception if record is not read in full
            Maybe useful for testing.
        """
        super().__init__(byteorder = byteorder)

        self.tolerant = tolerant
        self.open(filename,
                  extension = extension)

    def open(self,
             filename = None,
             extension = None,
             ext_exclude = ['xxx']):
        """
        Open the file.

        """

        if filename is not None:
            self.filename = os.path.expandvars(os.path.expanduser(filename))
        else:
            self.filename = ''

        compression_extensions = [''] + [os.path.extsep + ext for ext in ['gz', 'bz2', 'xz']]

        if extension is not None:
            if extension.startswith(os.path.extsep):
                extensions = [extension]
            else:
                extensions = [extension, os.path.extsep + extension]
            for cext,ext in itertools.product(compression_extensions, extensions):
                pattern = os.path.join(self.filename, '*' + ext)
                fp =  pattern + cext
                fx = glob.glob(fp)
                for fn in fx:
                    for x in ext_exclude:
                        if fn.startswith(x):
                            fn = ''
                            break
                    if os.path.isfile(fn):
                        break
                else:
                    fn = ''
                if os.path.isfile(fn):
                    self.filename = fn
                    break
            else:
                raise IOError("File not found.")
        else:
            for cext in compression_extensions:
                fp =  self.filename + cext
                fx = glob.glob(fp)
                for fn in fx:
                    if os.path.isfile(fn):
                        break
                else:
                    fn = ''
                if os.path.isfile(fn):
                    self.filename = fn
                    break
            else:
                raise IOError("File not found.")
        ext =  os.path.splitext(self.filename)[-1].strip(os.path.extsep)
        if ext == 'gz':
            self.compressed = True
            self.compress_mode = 'gz'
            self.file = gzip.GzipFile(self.filename, 'rb')
            pos = self.file.myfileobj.tell()
            self.file.myfileobj.seek(-4, os.SEEK_END)
            self.filesize = np.ndarray(1, dtype = "<u4",
                                       buffer = self.file.myfileobj.read(4))[0]
            self.file.myfileobj.seek(pos, os.SEEK_SET)
            self.stat = os.fstat(self.file.fileno())
        elif ext == 'bz2':
            self.compressed = True
            self.compress_mode = 'bz2'
            self.file = bz2.BZ2File(self.filename, 'rb', 2**24)
            # the following is lekely slow
            pos = self.file.tell()
            self.file.seek(0, os.SEEK_END)
            self.filesize = self.file.tell()
            self.file.seek(pos, os.SEEK_SET)
            self.stat = os.fstat(self.file.fileno())
        elif ext == 'xz':
            self.compressed = True
            self.compress_mode = 'xz'
            self.filesize = xz_file_size(self.filename)
            self.file = lzma.LZMAFile(self.filename, 'rb')
            self.stat = os.fstat(self.file.fileno())
        else:
            self.file = open(self.filename,'rb',-1)
            self.stat = os.fstat(self.file.fileno())
            self.filesize = self.stat.st_size
            self.compressed = False
        self._init()

    def _init(self):
        """
        Initialize the file position and data to empty.
        """
        self.fpos = 0
        self.pos = 0
        self.data = b''
        self.reclen = 0

    def eof(self):
        """
        Return 'EOF' status.

        True  if at or past end of data of last record.
        False otherwise.
        """
        return self.fpos  >= self.filesize

    def close(self):
        """Close the file."""
        self.file.close()

    def rewind(self):
        """Rewind the file."""
        self.file.seek(0,os.SEEK_SET)
        self._init()

    def load(self, *args, **kwargs):
        """Read a data record from file."""
        self._load(*args, **kwargs)

    def _load(self, size = None, offset = None):
        """Read in data of a record."""
        f = self.file

        if offset is not None:
            pos = offset
            f.seek(pos, os.SEEK_CUR)
        else:
            pos = 0

        if size is None:
            self.data = f.read()
            self.reclen = len(self.data)
        else:
            self.data = f.read(size)
            self.reclen = size

        if size is not None:
            if (not self.tolerant) and (self.reclen < size):
                raise FileReadError(self.filename, "Could not read requested number of bytes.")

        self.fpos += reclen
        self.pos = 0

    def assert_eof(self):
        """Throw exception if current position is not end of file.

        This can be use to deterime whether all data has been read,
        i.e., as a consistency check of the read data sizes.

        If the file is compressed, the initial file size cannot be
        used to determine the the size of the file."""
        if not self.eof() and not self.compressed:
            raise self.DataFileError(self.filename,
                                     self.filesize,
                                     self.fpos)

    # context manager interface

    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is None:
            self.close()
        return False

#=======================================================================

class FortranReader(DataReader, FortranSpecs):
    """
    Class for reading 'unformatted' Fortran binary files.

    This is in part based on a code from Laurens Keek (2010).

    Compressed files with .gz will be automatically uncompressed.
    This may fail, however, if the file is bigger tha 2GB or 4GB.

    Compressed files with .bz2 will be automatically uncompressed.
    This is currently rather inefficient because to determine the
    file size the entire stream need to be read first.

    Compressed files with .xz exension will be automatically
    uncompressed.  Uses a reasonably efficient code to determine total
    uncompressed file size, but does need to read end of compressed
    file on opening (could be slow on network-mounted or archived
    files).

    """

    def __init__(self,
                 *args,
                 byteorder = None,
                 reclen = 4,
                 verbose = False,
                 **kwargs):
        self._set_reclen(reclen)
        super().__init__(*args,
                         byteorder = byteorder,
                         **kwargs)
        if byteorder is None:
            if not self._check_byteorder():
                self._reverse_byteorder()
        self.verbose = verbose

    def _check_byteorder(self):
        """
        deterimine if file is opened in right endian
        """
        ok = True
        ok &= self.filesize >= 8
        pos0 = self.file.tell()
        if pos0 != 0:
            self.file.seek(0, os.SEEK_SET)
        reclen = self._read_reclen()
        ok &= (reclen >= 0) and (reclen <= self.filesize - 8)
        if ok:
            self.file.seek(reclen, os.SEEK_CUR)
            ok &= reclen == self._read_reclen()
        if pos0 != self.file.tell():
            self.file.seek(pos0, os.SEEK_SET)
        return ok

    def _set_byteorder(self, *args, **kwargs):
        super()._set_byteorder(*args, **kwargs)
        self.reclen_dtype = self.reclen_dtype.newbyteorder(self.byteorder)

    def _init(self):
        """
        Initialize the file position and data to empty.
        """
        super()._init()
        self.rpos = 0

    def iterrecords(self):
        """
        Return iterator over records in file in form of
        a readable recored.

        QUESTION - Do we need to copy the buffer?
        """

        class RecordReader(self.__class__):
            """
            class to read from buffer
            (one record)
            [could be more?]
            """

            def __init__(self,
                         data,
                         reclen,
                         pos,
                         byteorder = DataReader.default_byte_order,
                         *args,
                         **kwargs):

                self.data = data
                self.reclen = reclen
                self.pos = pos

                self._set_byteorder(byteorder = byteorder)

            def load(self):
                pass

        while not self.eof():
            self.load()
            yield RecordReader(
                self.data,
                self.reclen,
                self.pos,
                byteorder = self.byteorder)

    def _read_reclen(self):
        """Return the record length."""
        # fast python 2 code:
        # self.i4.data[0:4] = self.file.read(4)
        # return int(self.i4[0])

        # python 3 is more annoying
        # in particular, data now has chunks of size of the data type
        # and is a memory view rather than 'buffer'
        # the following seemed to have worked, though
        # b4 = np.ndarray(4, buffer=i4, dtype=np.byte) [add as class attribute]
        # b4.data[:]= memoryview(self.file.read(4))[:]
        # return int(b4)

        # the following works in both P2/3 but is slower.
        # it does create extra objects...
        return int(np.frombuffer(self.file.read(self.fortran_reclen),
                                 self.reclen_dtype))

    def skip(self, n = 1):
        """Read past n records (default: 1) without unpacking the data."""
        for i in range(n):
            self._load(False)

    def backspace(self, n = 1):
        """Read backward n records - will not work for compressed files."""
        f = self.file
        for i in range(n):
            f.seek(-self.fortran_reclen, os.SEEK_CUR)
            reclen = self._read_reclen()
            size = reclen + 2 * self.fortran_reclen
            f.seek(-size, os.SEEK_CUR)
            check = self._read_reclen()
            if (not self.tolerant) and (check != reclen):
                raise RecLenError(self.filename, "Header lenght does not match trailer length.")
            f.seek(-self.fortran_reclen, os.SEEK_CUR)
            self.fpos -= size
        self.rpos -= n
        self.pos = 0

    def _load(self, load = True, size = None, offset = None):
        """Read in data of a record or skip, and advance to next."""
        f = self.file
        if self.eof():
            raise EOFError()
        reclen = self._read_reclen()
        if reclen < 0:
            raise RecLenError(self.filename, "Negative record length.")

        if reclen == 0 and self.verbose:
            print(' [FortranRead] Warning: reclen == 0 in record {}'.format(self.rpos))

        if load:
            if offset is not None:
                if offset >= 0:
                    pos = offset
                else:
                    pos = reclen + offset
                f.seek(pos, os.SEEK_CUR)
            else:
                pos = 0
            if size is None:
                size = reclen
            else:
                if size >= 0:
                    size = min(size, reclen - pos)
                else:
                    size = max(reclen - pos + size, 0)
            self.data = f.read(size)
            self.reclen = size
            if size > len(self.data):
                raise ReadError(
                    '[_load] Could not read the requested number of bytes ({} vs {}).'.format(
                    len(self.data), size))
            remainder = reclen - pos - size
            if remainder > 0:
                f.seek(remainder, os.SEEK_CUR)
        else:
            f.seek(reclen, os.SEEK_CUR)
            self.data = b''
            self.reclen = 0

        check = self._read_reclen()
        if (not self.tolerant) and (check != reclen):
            raise RecLenError(self.filename, "Header lenght does not match trailer length.")

        self.fpos += reclen + 2 * self.fortran_reclen
        self.rpos += 1
        self.pos = 0

    def seek_noncorrupt(self, skip_empty = True, load = False):
        """Find next non-corrupt record candidate.

        This may not work for compressed files.

        This routine is still very preliminary."""
        # TODO - check next record and EOR?
        f = self.file
        tolerant, self.tolerant = self.tolerant, False
        while True:
            fpos = self.fpos
            try:
                f.seek(self.fpos, os.SEEK_SET)
                self._load(load)
            except:
                self.fpos += 1
            else:
                if skip_empty and self.fpos - fpos == 2 * self.fortran_reclen:
                    self.rpos -= 1
                    self.fpos = fpos + 1
                    continue
                break
        self.tolerant = tolerant
        if load is False:
            self.fpos = fpos
            f.seek(self.fpos, os.SEEK_SET)
