import re
import builtins

import numpy as np

# numpy / python type mape needed
_type_map = {'uint': 'int'}

# numpy types
_np_types = ('i1', 'i2', 'i4', 'i8',
          'u1', 'u2', 'u4', 'u8',
          'f2', 'f4', 'f8', 'f16',
          'c8', 'c16', 'c32',
          )

# function help generate new methods
def _set_method(
        cls = None,
        t = None,
        fn = None,
        name = None,
        parent = None,
        doc = None,
        extra_kw = None):
    if extra_kw is None:
        extra_kw = {}
    dt = np.dtype(t)
    dn = np.dtype(dt).name
    pt = re.findall('[a-z]+', dn)[0]
    pt = _type_map.get(pt, pt)
    f = builtins.__dict__[pt]
    if parent is not None:
        p = cls.__dict__[parent.format(**locals())]
    else:
        p = None
    _f = fn(f = f, p = p, dt = dt, **extra_kw)
    fn = name.format(**locals())
    _f.__qualname__ = "{}.{}".format(cls.__qualname__, fn)
    _f.__name__ = fn
    _f.__doc__ = doc.format(**locals())
    setattr(cls, fn, _f)
