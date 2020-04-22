"""
Provide various utilties not found in Python 2.7 or numpy 1.5
"""

import functools
import operator
import copy
import datetime
import os
import struct
import time
import logging
import sys
from collections import Iterable
from functools import partial
from glob import glob
import numpy as np

def cumsum(list,seed=0):
    """Return List with cumulative sum starting at seed (exclusive)."""
    return [x for x in cumsum_iter(list,seed)]

def cumsum0(list,seed=0):
    """Return List with cumulative sum starting at seed (inclusive)."""
    return [x for x in cumsum0_iter(list,seed)]

def cumsum_iter(list,seed=0):
    """Cumulative sum starting at seed (default 0)."""
    t = seed
    for l in list:
        t += l
        yield t

def cumsum0_iter(list,seed=0):
    """Iterator for cumulative sum starting at seed (default 0)."""
    t = seed
    yield t
    for l in list:
        t += l
        yield t

def cumsum0_range_iter(list,seed=0):
    """Iterator for cumulative sum duets for slicing starting at seed (default 0)."""
    t1 = seed
    for l in list:
        t0 = t1
        t1 += l
        yield t0,t1

def cumsum0_enum_range_iter(list,seed=0):
    """Iterator for cumulative sum and counter triples for slicing starting at seed (default 0)."""
    t1 = seed
    for i,l in enumerate(list):
        t0 = t1
        t1 += l
        yield i,t0,t1

def prod(seq):
    """Product of a sequence."""
    if not np.iterable(seq):
        seq = (seq,)
    return functools.reduce(operator.mul, seq, 1)


def contract(a,
             sequence,
             axis=0,
             dimension = None):
    """
    Contract array from list.

    typical use would be:
        categories, members = np.unique(X.index, return_inverse = True)
        result = contract(X.values, members)

    replacing:
        result = np.zeros((..., len(categories), ...), dtype = np.float64)
        values = X.values
        for i,j in enumerate(members):
             result[..., j,...] += values[..., i,...]
    """
    shape = np.array(a.shape)
    ii = [slice(i) for i in shape]
    jj = copy.deepcopy(ii)
    if axis == -1:
        axis += a.ndim
    axis_dimension = np.amax(sequence) + 1
    if dimension is None:
        dimension = axis_dimension
    else:
        assert dimension >= axis_dimension, "Target dimension too small."
    shape[axis] = dimension
    out = np.zeros(shape)
    for i,j in enumerate(sequence):
        jj[axis] = j
        ii[axis] = i
        out[jj] += a[ii]
    return out

def project(a,
            values,
            axis=0,
            return_values = False,
            return_inverse = False):
    """
    Project array onto values.
    """
    k, kinv = np.unique(values, return_inverse = True)
    p = [contract(a, kinv, axis = axis)]
    if return_values:
        p.append(k)
    if return_inverse:
        p.append(kinv)
    if len(p) == 1:
        p = p[0]
    else:
        p = tuple(p)
    return p


def isinslice(index, Slice):
    """
    Determine whether index is part of a slice.
    """
    start, stop, stride = Slice.indices(index + 1)
    if (index - start) % stride != 0:
        return False
    if stride < 0:
        return (start >= index > stop)
    else:
        return (start <= index < stop)

def bool2sign(true):
    return 1 if true else -1

def sign(integer):
    return (integer > 0) - (integer < 0)

def bool2slice(true):
    return slice(None,None,1) if true else slice(None,None,-1)

class UTCFormatter(logging.Formatter):
    converter = time.gmtime
    def format(self, record):
        record.nameb = '[{:s}]'.format(record.name)
        return super(UTCFormatter, self).format(record)

class Slice(object):
    """
    Slice iterator object.
    """
    def __init__(self, *args, **kwargs):
        """
        Construct from slice indices or slice object.  Provide
        optional object size.
        """
        if len(args) == 1 and isinstance(args[0], slice):
            self.slice = args[0]
        else:
            self.slice = slice(*args)
        self.size = kwargs.pop('size', None)
        assert len(kwargs) == 0
    def __iter__(self):
        """
        Slice object iterator.
        """
        if self.size is None:
            self.size = max(self.slice.start, self.slice.stop) + 1
        xslice = self.slice.indices(self.size)
        for i in range(*xslice):
            yield i
    def iter(self, size = None):
        """
        Return iterator with defined object size.
        """
        size = self.size
        for i in self.__iter__():
            yield i
        self.size = size


class ClearCache(object):

    def clear_cache(self, debug = False):
        """
        delete cached variable
        """
        cls = type(self)
        try:
            cls._del_attr
        except:
            import types
            import inspect
            my_name = inspect.stack()[0][3]
            cls._del_attr = []
            cls._del_meth = []
            for d in dir(cls):
                if d == my_name:
                    continue
                D = getattr(cls, d)
                if isinstance(D, CachedAttribute):
                    cls._del_attr += [d]
                if isinstance(D, types.FunctionType) and hasattr(D, 'clear'):
                    cls._del_meth += [d]

        for d in cls._del_attr:
            try:
                del self.__dict__[d]
            except KeyError:
                if debug:
                    print(' [DEBUG] "{}" not found'.format(d))
                pass

        for d in cls._del_meth:
            getattr(cls, d).__call__(self, clear = True)


def cachedmethod(method, *args, **kwargs):
    """
    Decorator to compute a quantity on first call only.

    The data is stored in '_'+method.__name__

    Use:
    class X(parent_object):
        ...
        @cachedmethod
        def method(parameters):
            ...

    Use function clear_cachedmethod to delete data field.
    if you have an instance x of your class X:

    >>> x = X(init_parameters)

    you would call

    >>> clear_cachedmethod(x.method)

    You can also call

    >>> x.X(clear = True)

    To force recalculation, call

    >>> x.X(recalc = True)

    """
    key = '_' + method.__name__
    args_name = key + '_args'
    kwargs_name = key + '_kwargs'

    def cached_method(self, *args, **kwargs):
        """
        Method only to be called first time or when parameters change.
        """
        clear_data = kwargs.pop('clear', False)
        if clear_data:
            clear(self)
            return
        recalc = kwargs.pop('recalc', False)
        if recalc:
            clear(self)
        d = self.__dict__
        reload = key not in d
        if ((method.__code__.co_argcount > 1) or
            (method.__code__.co_flags & 0x0c > 0)):
            if not reload:
                if d.get(args_name, None) != args:
                    d[args_name] = args
                    reload = True
                if d.get(kwargs_name, None) != kwargs:
                    d[kwargs_name] = kwargs
                    reload = True
            else:
                if not args_name in d:
                    d[args_name] = args
                if not kwargs_name in d:
                    d[kwargs_name] = kwargs
        if reload:
            d[key] = method(self, *args, **kwargs)
        return d[key]

    def clear(self):
        """
        Clear storage.  Requires class instance passed explicitly.

        >>> x.method.clear(x)
        """
        d = self.__dict__
        d.pop(key, None)
        d.pop(args_name, None)
        d.pop(kwargs_name, None)

    cached_method.__dict__.update(method.__dict__)
    cached_method.__dict__['clear'] = clear
    cached_method.__dict__['method'] = method.__name__
    if method.__doc__ is not None:
        cached_method.__doc__ = method.__doc__ + '\n' + cached_method.__doc__
    cached_method.__name__ = method.__name__
    cached_method.__module__ = getattr(method, '__module__')
    return cached_method

def clear_cachedmethod(method):
    """
    Clear the stored data for a method created with the @firstcall decorator.

    >>> clear_cachedmethod(x.method)
    """
    method.clear(method.__self__)


class CachedAttribute(object):
     """
    Computes attribute value and caches it in the instance.

    Source:
    http://stackoverflow.com/questions/3237678/how-to-create-decorator-for-lazy-initialization-of-a-property

    Reference as given in source:
    http://code.activestate.com/recipes/276643-caching-and-aliasing-with-descriptors/

    Use 'del inst.myMethod' to clear cache.

    Note that if this depends on other cached attribute, those will
    have to be clear indendently and directly.
    """
     def __init__(self, method, name = None):
         self.method = method
         self.name = name or method.__name__
     def __get__(self, obj, objtype):
         if obj is None:
             return self
         elif self.name in obj.__dict__:
             return obj.__dict__[self.name]
         else:
             value = self.method(obj)
             obj.__dict__[self.name] = value
             return value
     def __set__(self, obj, value):
         raise AttributeError("Cannot assign to " + self.name + '.')
     def __delete__(self, obj):
         try:
             del obj.__dict__[self.name]
         except KeyError:
             pass

class CachedClassAttribute(object):
    """
    Computes class attribute value and caches it in the class.

    Source:
    http://stackoverflow.com/questions/3237678/how-to-create-decorator-for-lazy-initialization-of-a-property

    Reference as given in source:
    http://code.activestate.com/recipes/276643-caching-and-aliasing-with-descriptors/

    Use 'del MyClass.myMethod' to clear cache.
    """
    # has not been tested ...
    # seems to conflict with autoreload of IPython
    def __init__(self, method, name = None):
        self.method = method
        self.name = name or method.__name__
    def __get__(self, obj, objtype):
        if self.name in objtype.__dict__:
            return getattr(objtype, self.name)
        else:
            value = self.method(objtype)
            setattr(objtype, self.name, value)
            return value
    def __set__(self, obj, value):
        raise AttributeError("Cannot assign to " + self.name + '.')
    def __delete__(self, objtype):
        # does this really pass type?
        # del objtype.__dict__[self.name]
        # maybe the following works as well
        delattr(objtype, self.name)


def Property(func):
    """
    Use:
    class Person(object):
    @Property
    def name():
        doc = "The person name"

        def fget(self):
            return _name

        def fset(self, name):
            self._name = name

        def fdel(self):
            del self.last_name

        return locals()
    """
    return property(**func())

from types import MethodType
def make_cached_attribute(self, func, name = None, doc = None, args = None, kw = None):
    """
    Add cached attribute to class dynamaically.

    EXAMPLE:
        (call in __init__)

        def var(self, idx):
            return np.array([x.output[idx] for x in self.data])

        make_cached_attribute(self.__class__,
                              functools.partial(var,idx=21),
                              'xh','central XH')
    """
    if args is None:
        args = list()
    if kw is None:
        kw = dict()
    setattr(self, name + '_kw', kw)
    setattr(self, name + '_args', args)
    def f(self):
        kw = self.__getattribute__(name + '_kw')
        args = self.__getattribute__(name + '_args')
        return func(self, *args, **kw)
    f.__doc__  = doc
    f.__self__.__class__ = self.__class__
    f.__func__  = f
    f.__self__  = None
    f = CachedAttribute(f, name)
    # def __get__(self, instance, owner):
    #     return MethodType(self, instance, owner)
    # f.__get__ = __get__
    setattr(self.__class__, name, f)

class OutFile(object):
    """
    Contex Manager: Open file if filename is given or use file.
    """
    def __init__(self,
                 outfile = None,
                 silent = False,
                 overwrite = False):
        """
        open `stdout` if file does not exist.
        """
        # self.setup_logger(silent = silent)
        if outfile is None:
            self.f = sys.stdout
            self.open = False
            return
        self.open = isinstance(outfile, file)
        if not self.open:
            filename = os.path.expanduser(os.path.expandvars(outfile))
            assert overwrite or not os.path.exists(filename)
            f = open(filename,'w')
        else:
            f = outfile
        self.f = f

    def __enter__(self):
        return self.f

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.open:
            self.f.close()
        # self.close_logger()

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
        f.seek(-8, os.SEEK_END)
        bkwd_size = np.ndarray((), dtype = "<u4", buffer = f.read(4))
        # print(bkwd_size)

        # 12 byte footer (4-CRC,4-size,2-flags,2-YZ)
        # 4 * (backward_size + 1) is start of index
        # index starts with 0x00 flag, last 4 byte arce CRC
        f.seek(-12 - 4 * (bkwd_size + 1), os.SEEK_END)
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


class MultiLoop(object):
    """
    Provides multi_loop method to loop over all iterable parameters
    except strings.

    Use:
        class X( ..., MultiLoop, ...)

    Call:
        self.multi_loop(self.method_to_run, *args, **kwargs)

    Parameters:
        loop_descend:
            a keyword parameter that decides what to do with nested
            iterables
        no_resolve:
            overwrite classes not to resolve
        no_resolve_add
            add classes not to resolve

    LIMITATIONS:
        In contrast, dict[1] in second level will be preserved.
        (largely, because it cannot be resolved in a reasonable way in
        the first place)

    TODO:
        resove nested lists in lists in parameter items
          {mass: (1,2)} --> {mass: 1}, {mass: 2}
          {('A','B'): X} --> {'A': X}, {'B': X}
        maybe should have separate keywords for
          list/tuple and dictionary resolution (depth)
    """
    no_resolve_ = (str, set)
    class _multi_loop_container(object):
        def __init__(self, item, level = 0):
            self.item = item
            self.level = level
    def multi_loop(self, method, *args, **kwargs):
        """
        Loops over all iterable parameters except strings and sets.
        """
        kwargs = kwargs.copy()
        descend = kwargs.pop('loop_descend', 1)
        no_resolve = kwargs.pop('no_resolve', self.no_resolve_)
        no_resolve += kwargs.pop('no_resolve_add', tuple())
        kwargs_new = kwargs.copy()
        args_new = list(args)

        if descend <= 0:
            return [method(*args_new, **kwargs_new)]
        kwargs_new['loop_descend'] = descend
        kwargs_new['no_resolve'] = no_resolve
        result = []
        for iv,v in enumerate(args):
            if isinstance(v, no_resolve):
                continue
            level = 1
            if isinstance(v, self._multi_loop_container):
                if v.level == descend:
                    continue
                level = v.level + 1
                v = v.item
            if isinstance(v, dict):
                if len(v) <= 1:
                    continue
                for k,i in v.items():
                    args_new[iv] = {k:i}
                    result += self.multi_loop(method, *args_new, **kwargs_new)
                return result
            if isinstance(v, Iterable):
                for i in v:
                    args_new[iv] = self._multi_loop_container(i, level)
                    result += self.multi_loop(method, *args_new, **kwargs_new)
                return result
        for kw,v in kwargs.items():
            if isinstance(v, no_resolve):
                continue
            level = 1
            if isinstance(v, self._multi_loop_container):
                if v.level == descend:
                    continue
                level = v.level + 1
                v = v.item
            if isinstance(v, dict):
                if len(v) <= 1:
                    continue
                for k,i in v.items():
                    kwargs_new[kw] = {k:i}
                    result += self.multi_loop(method, *args_new, **kwargs_new)
                return result
            if isinstance(v, Iterable):
                for i in v:
                    kwargs_new[kw] = self._multi_loop_container(i,level)
                    result += self.multi_loop(method, *args_new, **kwargs_new)
                return result
        # get rid of containers
        for iv,v in enumerate(args_new):
            if isinstance(v, self._multi_loop_container):
                args_new[iv] = v.item
        for kw,v in kwargs_new.items():
            if isinstance(v, self._multi_loop_container):
                kwargs_new[kw] = v.item
        return [method(*args_new, **kwargs_new)]
    @staticmethod
    def clean(kwargs, extra = None):
        """
        clean out MultiLoop kw arguments
        """
        kw = kwargs.copy()
        if extra is not None:
            if isinstance(extra, str):
                extra = (extra,)
        else:
            extra = tuple()
        extra += ('loop_descend', 'no_resolve', 'no_resolve_add')
        for x in extra:
            kw.pop(x, None)
        return kw

def loopmethod(descend = 1,
               no_resolve = (str, set),
               no_resolve_add = tuple(),
               ):
    """
        Decorator to compute a looped method.

        Use:
        @loopmethod(kwargs)
        method_to_loop

        If descend is False, stope at first level, otherwise descend
        down nested lists, sets, and tuples.

        Call:
        self.method_to_loop(*args,**kwargs)

        Returns list of results.

        kwargs:
        descend:
            a keyword parameter that decides what to do with nested
            iterables
        no_resolve:
            overwrite classes not to resolve
        no_resolve_add
            add classes not to resolve
    """
    if descend <= 0:
        return [method(self, *args, **kwargs)]
    if no_resove is None:
        no_resolve = tuple
    if not isinstance(no_resolve, Iterable):
        no_resolve = tuple((no_resolve,))
    if not isinstance(no_resolve_add, Iterable):
        no_resolve_add = tuple((no_resolve_add,))
    no_resolve = tuple(no_resolve) + tuple(no_resolve_add)

    def loop_method(method):
        """
        Decorator to compute a looped method.

        Use:
        @loopedmethod
        method_to_loop

        Call:
        self.method_to_loop(*args,**kwargs)

        """
        class _container(object):
            def __init__(self, item, level):
                self.item = item
                self.level = level
        # @wraps(method)
        def looped_method(self, *args, **kwargs):
            """
            Loop over all Iterables in *args and **kwargs except strings and sets.
            """
            kwargs_new = kwargs.copy()
            args_new = args.copy()
            result = []
            for iv,v in enumerate(args):
                if isinstance(v, no_resolve):
                    continue
                level = 1
                if isinstance(v, _container):
                    if v.level == descend:
                        continue
                    level = v.level + 1
                    v = v.item
                if isinstance(v, dict):
                    if len(v) <= 1:
                        continue
                    for k,i in v.items():
                        args_new[iv] = {k:i}
                        result += looped_method(self, *args_new, **kwargs_new)
                    return result
                if isinstance(v, Iterable):
                    for i in v:
                        args_new[iv] = _container(i, level)
                        result += looped_method(self, *args_new, **kwargs_new)
                    return result
            for kw,v in kwargs.items():
                if isinstance(v, no_resolve):
                    continue
                level = 1
                if isinstance(v, _container):
                    if v.level == descend:
                        continue
                    level = v.level + 1
                    v = v.item
                if isinstance(v, dict):
                    if len(v) <= 1:
                        continue
                    for k,i in v.items():
                        kwargs_new[kw] = {k:i}
                        result += looped_method(self, *args_new, **kwargs_new)
                    return result
                if isinstance(v, Iterable):
                    for i in v:
                        kwargs_new[kw] = _container(i, level)
                        result += looped_method(self, *args_new, **kwargs_new)
                    return result
            # get rid of containers
            for iv,v in enumerate(args_new):
                if isinstance(v, _container):
                    args_new[iv] = v.item
            for kw,v in kwargs_new.items():
                if isinstance(v, _container):
                    kwargs_new[kw] = v.item
            return [method(self, *args_new, **kwargs_new)]
        looped_method.__dict__.update(method.__dict__)
        looped_method.__dict__['method'] = method.__name__
        if method.__doc__ is not None:
            looped_method.__doc__ = method.__doc__ + '\n' + looped_method.__doc__
        looped_method.__name__ = method.__name__
        looped_method.__module__ = getattr(method, '__module__')
        return looped_method
    return loop_method

# class test(object):
#     def __init__(self, *args, **kwargs):
#         X = self.x(*args, **kwargs)
#         print(X)
#     @loopmethod(3)
#     def x(self,*args, **kwargs):
#         print(args,kwargs)

class test(MultiLoop):
    def __init__(self, *args, **kwargs):
        X = self.multi_loop(self.x, *args, loop_descend = 2, **kwargs)
        print(X)
    def x(self,*args, **kwargs):
        base = kwargs.setdefault('base',None)
        print(list(base.values()))
        print(args,kwargs)


def float2str(f, precision = 13):
    """
    Use g format but add '.' to be compatible with KEPLER
    """
    s = ("{:."+str(precision)+"g}").format(f)
    if ((s.find('.') == -1) and
        (s.find('e') == -1)):
        s += '.'
    if ((s.find('.') == -1) and
        (s.find('e') != -1)):
        s = s.replace('e','.e')
    return s


def bit_count(x):
    assert x>=0, 'negative numbers have infinite number of leading bits'
    count = 0
    bit = 1
    while bit <= x:
        count += int((x & bit) > 0)
        bit <<= 1
    return count

def queue_processor(input, output, params):
    """
    Worker thread to process data from queue.

    Assume input is multiprocessing.JoinableQueue or Queue.Queue.

    call signature
      queue_processor(input, output, params)

    input
      is a Queue provides get() and task_done()

    output
      is a Queue that provides put()

    params should be a dictionary that contains
      data
        basic initialization data
        could be a large data set to operate on
      processor
        a class that is initialzed with data
        __init__(data)
        it is called to provide results that are put in out queue
        __call__(task)
    """
    # do we really need a dictionary?
    processor = params.get('processor')(
        params.get('data', None))
    # just to make sure we remove unnecessay references
    # as we may have many copies in parallel processes
    del params
    for task in iter(input.get, 'STOP'):
        output.put(processor(task))
        input.task_done()
    input.task_done()


def stuple(*args):
    """
    convert string to tuple with one element,
    list to tuple, None to empty tuple,
    leave tuple unchanged.
    """
    out = ()
    for a in args:
        if isinstance(a, str):
            s = (' ' + a + ' ').splitlines()
            s[0] = s[0][1:]
            s[-1] = s[-1][:-1]
            out += tuple(s)
        elif a is not None:
            assert isinstance(a, (list, tuple, np.ndarray))
            for b in a:
                out += stuple(b)
    return out


import physconst
def ergs2mbol(ergs):
    return +4.77 - log10(ergs / physconst.XLSUN) * 2.5

# the following now is a singleton metaclass
class SMeta(type):
    """
    Usage:
       class X(Y, metaclass = MetaSingletonHash)
    """
    def __call__(*args):
        cls = args[0]
        key = args
        try:
            cache = cls._cache
        except:
            cache = dict()
            cls._cache = cache
        try:
            obj = cache[key]
        except:
            obj = type.__call__(*args)
            cache[key] = obj
        return obj

# the following now is a singleton metaclass
class MetaSingletonHash(type):
    """
    Singleton metaclass based on hash

    First creates object to be able to test hash.

    If same hash is found, return old object and discard new one,
    otherwise return old one.

    Usage:
       class X(Y, metaclass = MetaSingletonHash)

    class X needs to provide a __hash__ function
    """
    def __call__(*args, **kwargs):
        cls = args[0]
        try:
            cache = cls._cache
        except:
            cache = dict()
            cls._cache = cache
        obj = type.__call__(*args, **kwargs)
        key = (cls.__name__, obj.__hash__())
        return cache.setdefault(key, obj)

import gzip
import bz2
import lzma

def text_file(filename = None,
              mode = None,
              compress = None,
              return_filename = False,
              return_compress = False,
          ):
    if mode is None:
        if compress:
            mode = 'w'
        else:
            mode = 'r'
    if filename:
        filename = os.path.expandvars(os.path.expanduser(filename))
        if compress == True:
            compress = 'xz'
        if compress:
            if not (filename.endswith('.gz') or
                    filename.endswith('.bz2') or
                    filename.endswith('.xz')):
                if not compress.startswith('.'):
                    compress = '.' + compress
                filename += compress
        else:
            compress = False
        if filename.endswith('.gz'):
            fout = gzip.open(filename, mode + 't', encoding = 'ASCII')
        elif filename.endswith('.bz2'):
            fout = bz2.open(filename, mode + 't', encoding = 'ASCII')
        elif filename.endswith('.xz'):
            fout = lzma.open(filename, mode + 't', encoding = 'ASCII')
        else:
            fout = open(filename, mode + 't', -1)
    else:
        fout = sys.stdout
        compress = False
        filename = sys.stdout.name
    out = [fout]
    if return_filename:
        out += [filename]
    if return_compress:
        out += [compress]
    if len(out) == 1:
        return out[0]
    return tuple(out)

def close_text_file(fout):
    if fout != sys.stdout:
        fout.close()

class TextFile(object):
    def __init__(self,
                 *args,
                 return_filename = False,
                 return_compress = False,
                 **kwargs):
        self.file, self.filename, self.compress = text_file(
            *args,
            return_filename = True,
            return_compress = True,
            **kwargs)
        self.return_filename = return_filename
        self.return_compress = return_compress
    def write(self, *args, **kwargs):
        self.file.write(*args, **kwargs)
    def writable(self):
        return self.file.writable()
    def writelines(self, *args, **kwargs):
        self.file.writelines(*args, **kwargs)
    def read(self, *args, **kwargs):
        return self.file.read(*args, **kwargs)
    def readable(self):
        return self.file.readable()
    def readline(self, *args, **kwargs):
        return self.file.readline(*args, **kwargs)
    def readlines(self, *args, **kwargs):
        return self.file.readlines(*args, **kwargs)
    def tell(self):
        return self.file.tell()
    def seekable(self, *args, **kwargs):
        return self.file.seekable(*args, **kwargs)
    def seek(self, *args, **kwargs):
        self.file.seek(*args, **kwargs)
    @property
    def encoding(self):
        return self.file.encoding
    def close(self):
        if self.file != sys.stdout:
            self.file.close()
    def __enter__(self):
        out = [self.file]
        if self.return_filename:
            out += [self.filename]
        if self.return_compress:
            out += [self.compress]
        if len(out) == 1:
            return out[0]
        return tuple(out)
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

def iterable(x):
    """
    convert things to an iterable, but omit strings

    May need to add other types.
    """
    if isinstance(x, str):
        x = (x,)
    if isinstance(x, np.ndarray) and len(x.shape) == 0:
        x = (x,)
    if not isinstance(x, (Iterable, np.ndarray)):
        x = (x,)
    return x

def np_array(x):
    if isinstance(x, np.ndarray) and len(x.shape) > 0:
        return x
    return np.array(iterable(x))

def is_iterable(x):
    """
    return whether is a true iterable incluidng numpy.ndarra, but not string
    """
    if isinstance(x, np.ndarray) and len(x.shape) == 0:
        return False
    return isinstance(x, (Iterable, np.ndarray)) and not isinstance(x, str)

import contextlib

@contextlib.contextmanager
def environ(env):
    assert isinstance(env, dict)
    save = dict()
    for key, val in env.items():
        save[key] = os.environ.get(key, None)
        os.environ[key] = val
    yield
    for key, val in save.items():
        if val is None:
            del os.environ[key]
        else:
            os.environ[key] = val

def walk_files(path, ignore = None):
    # todo - make ignore a RE match
    for dirpath, dirs, files in os.walk(path):
        for f in files:
            yield os.path.join(dirpath, f)
        for d in iterable(ignore):
            try:
                dirs.remove(d)
            except:
                pass


def touch(filename, verbose = False, timestamp = None):
    """
    Update file to current date or timestamp (seconds), create empty
    file if it does not exist.
    """
    xtime = None
    if timestamp is not None:
        xtime = (timestamp,)*2
    if verbose:
        print('Touching {}'.format(filename))
    try:
        os.utime(filename, xtime)
    except:
        open(filename, 'a').close()


@contextlib.contextmanager
def chdir(path = None):
    """
    Context managet to work in provided directory.
    """
    cwd = os.getcwd()
    try:
        if path is not None:
            os.chdir(path)
        yield
    finally:
        os.chdir(cwd)


# class NestedDict(dict):
#     """
#     >>> eggs = NestedDict()
#     >>> eggs[1][2][3][4][5]
#     {}
#     >>> eggs
#     {1: {2: {3: {4: {5: {}}}}}}
#     """
#     def __getitem__(self, key):
#         if key in self: return self.get(key)
#         return self.setdefault(key, NestedDict())

# class MyDict(dict):
#     def __missing__(self, key):
#         t = self[key] = MyDict()
#         return t

# MyDict = lambda: collections.defaultdict(MyDict)

def make_nan(payload = 0, sign = False, quiet = True):
    """IEEE 754-2008"""
    if (not quiet) and (payload == 0):
        raise AttributeError('NaN cannot be signaling with payload 0.')
    if not 0 <= payload <= 0x0007ffffffffffff:
        raise AttributeError('Payload outside range.')
    payload |= (0x7ff0000000000000 |
                (sign * 0x8000000000000000) |
                (quiet *  0x0008000000000000))
    return struct.unpack('d', struct.pack('Q', payload))[0]

def scan_nan(nan):
    """IEEE 754-2008"""
    nan = struct.unpack('Q', struct.pack('d', nan))[0]
    if ((nan & 0x7ff0000000000000 != 0x7ff0000000000000) or
        (nan & 0x000fffffffffffff == 0)):
        raise AttributeError('Attribute is not a NaN.')
    sign = nan & 0x8000000000000000 != 0
    quiet = nan & 0x0008000000000000 != 0
    payload = nan & 0x0007ffffffffffff
    return payload, sign, quiet
