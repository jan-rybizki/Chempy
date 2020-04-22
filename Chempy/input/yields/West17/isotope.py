"""
Collection of routines for isotopes, most prominently the Ion class.
"""

# contains all the isotope routines

import string
import re
import copy
import collections
import builtins

from collections.abc import Iterable, Mapping

import numpy as np

from utils import MetaSingletonHash

# # maybe this should be an enumeration???
# TODO - this class draft needs more work for tu
# class _elements(object):
#     _data = (
#     'nt',  'h',   'he',  'li',  'be',  'b',   'c',   'n',
#     'o',   'f',   'ne',  'na',  'mg',  'al',  'si',  'p',
#     's',   'cl',  'ar',  'k',   'ca',  'sc',  'ti',  'v',
#     'cr',  'mn',  'fe',  'co',  'ni',  'cu',  'zn',  'ga',
#     'ge',  'as',  'se',  'br',  'kr',  'rb',  'sr',  'y',
#     'zr',  'nb',  'mo',  'tc',  'ru',  'rh',  'pd',  'ag',
#     'cd',  'in',  'sn',  'sb',  'te',  'i',   'xe',  'cs',
#     'ba',  'la',  'ce',  'pr',  'nd',  'pm',  'sm',  'eu',
#     'gd',  'tb',  'dy',  'ho',  'er',  'tm',  'yb',  'lu',
#     'hf',  'ta',   'w',  're',  'os',  'ir',  'pt',  'au',
#     'hg',  'tl',  'pb',  'bi',  'po',  'at',  'rn',  'fr',
#     'ra',  'ac',  'th',  'pa',  'u',   'np',  'pu',  'am',
#     'cm',  'bk',  'cf',  'es',  'fm',  'md',  'no',  'lr',
#     'rf',  'db',  'sg',  'bh',  'hs',  'mt',  'ds',  'rg',
#     'cn',  'ut',  'fl',  'up',  'lv',  'us',  'uo')

#     def __getitem__(self, index):
#         if isinstance(index, str):
#             try:
#                 index = int(index)
#             except:
#                 pass
#         try:
#             return self._data[index]
#         except:
#             pass
#         try:
#             return self.index(index)
#         except:
#             raise IndexError()
#     def __len__(self):
#         return len(self._data)
#     def index(self, item):
#         return self._data.index(item.lower())
#     def __call__(self, value):
#         return self.__getitem__(value)
#     def __getattr__(self, attr):
#         try:
#             return self.index(attr)
#         except:
#             raise AttributeError()
# elements = _elements()

# should be Uus Uuo but I will not use 3 letters
# as this breaks too many things

# class _Elements(_elements):
#     _data = [e.capitalize() for e in _elements._data]
#     _data[0] = 'nt'
#     _data = tuple(_data)
#     def index(self, value):
#         try:
#             return self._data[value]
#         except ValueError:
#             return self._data[value.capitalize()]
# Elements = _Elements()

elements = (
    'nt',  'h', 'he', 'li', 'be',  'b',  'c', 'n',   'o',  'f',
    'ne', 'na', 'mg', 'al', 'si',  'p',  's', 'cl', 'ar',  'k',
    'ca', 'sc', 'ti',  'v', 'cr', 'mn', 'fe', 'co', 'ni', 'cu',
    'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr', 'rb', 'sr',  'y',
    'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', 'cd', 'in',
    'sn', 'sb', 'te',  'i', 'xe', 'cs', 'ba', 'la', 'ce', 'pr',
    'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy', 'ho', 'er', 'tm',
    'yb', 'lu', 'hf', 'ta',  'w', 're', 'os', 'ir', 'pt', 'au',
    'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn', 'fr', 'ra', 'ac',
    'th', 'pa', 'u',  'np', 'pu', 'am', 'cm', 'bk', 'cf', 'es',
    'fm', 'md', 'no', 'lr', 'rf', 'db', 'sg', 'bh', 'hs', 'mt',
    'ds', 'rg', 'cn', 'nh', 'fl', 'mc', 'lv', 'ts', 'og', 'ue')

Elements = [element.capitalize() for element in elements]
Elements[0] = 'nt'
Elements = tuple(Elements)

z2Name = (
    "Neutron", "Hydrogen", "Helium", "Lithium", "Beryllium",
    "Boron", "Carbon", "Nitrogen", "Oxygen", "Fluorine",
    "Neon", "Sodium", "Magnesium", "Aluminium", "Silicon",
    "Phosphorus", "Sulfur", "Chlorine", "Argon", "Potassium",
    "Calcium", "Scandium", "Titanium", "Vanadium", "Chromium",
    "Manganese", "Iron", "Cobalt", "Nickel", "Copper",
    "Zinc", "Gallium", "Germanium", "Arsenic", "Selenium",
    "Bromine", "Krypton", "Rubidium", "Strontium", "Yttrium",
    "Zirconium", "Niobium", "Molybdenum", "Technetium", "Ruthenium",
    "Rhodium", "Palladium", "Silver", "Cadmium", "Indium",
    "Tin", "Antimony", "Tellurium", "Iodine", "Xenon",
    "Caesium", "Barium", "Lanthanum", "Cerium", "Praseodymium",
    "Neodymium", "Promethium", "Samarium", "Europium", "Gadolinium",
    "Terbium", "Dysprosium", "Holmium", "Erbium", "Thulium",
    "Ytterbium", "Lutetium", "Hafnium", "Tantalum", "Tungsten",
    "Rhenium", "Osmium", "Iridium", "Platinum", "Gold",
    "Mercury", "Thallium", "Lead", "Bismuth", "Polonium",
    "Astatine", "Radon", "Francium", "Radium", "Actinium",
    "Thorium", "Protactinium", "Uranium", "Neptunium", "Plutonium",
    "Americium", "Curium", "Berkelium", "Californium", "Einsteinium",
    "Fermium", "Mendelevium", "Nobelium", "Lawrencium", "Rutherfordium",
    "Dubnium", "Seaborgium", "Bohrium", "Hassium", "Meitnerium",
    "Darmstadtium", "Roentgenium", "Copernicium", "Nihonium", "Flerovium",
    "Moscovium", "Livermorium", "Tennessine", "Oganesson", "(Ununennium)",
    )
z2name = tuple((name.lower() for name in z2Name))


# el2z provides easy conversion from element symbol to charge number
el2z = {el:i for i,el in enumerate(elements)}
el2z.update({el:i for i,el in enumerate(Elements[1:],start=1)})

def el2name(sym):
    return z2name[el2z[sym]]

class Ion(object, metaclass = MetaSingletonHash):
    """
    This is a generic class for 'ions'.
    It is designed to handel elements, isotioes, isomers, and isobars.
    Isomers have an energy level appended by 'e' or 'E';
        if the level is '0' then it should not be printed.
    It provides also an "index" field (idx) that can be used for sorting.
    __add__ and __sub__ should help with networks.
    R([[I],O]) allows to formulate reactions
        will return ground state nuclei
    Use (solely)
        A = ... for 'isobars'
        Z = ... for 'elements'
        N = ... for 'isotones'
    To specify isotopes you need to supply at least 2 of A, Z, N;
        if you specify three, they must fulfil A = N + Z
    To specify an isomer, you have provide additionally E = ...
        Use E = 0 to specify an ismer in the gs.

    g    = ( 0, 0, 0, 1)
    e-   = ( 7,-1, 0, 1)
    e+   = ( 7,+1, 0, 2)
    void = (-1, 0, 0, 0)

    We always put E=1, i.e., m1, as 'm'
    Default use 'mE' to indicate exited state E.
    We destinguish isotopes and isomers - Al26g is not Al26.
    Use EANY to specify generalised excited state, signified by '*'.
    """
    # in the future we want to use these
    # should flags be first or last?
    # maybe last?
    AMAX = 512
    ZMAX = 128 # -7 ... 120
    ZLIM = 120
    EMAX = 256
    FMAX = 128 # flags - 7 bit
    # leave the sign bit for int64

    FMUL = 1
    EMUL = FMUL * FMAX
    AMUL = EMUL * EMAX
    ZMUL = AMUL * AMAX
    BMUL = ZMUL * ZMAX
    # what lies Beyond ... (extra information for subclasses)
    # if these are used, the index no longer fits into int64

    # bit 0-2
    F_BOSON   = 0
    F_ELEMENT = 1
    F_ISOTOPE = 2
    F_ISOMER  = 3
    F_ISOBAR  = 4
    F_ISOTONE = 5
    F_HADRON  = 6 # currently not used
    F_LEPTON  = 7
    F_GROUP_MASK = 7

    # bits 3-6 remain for other implementation
    F_OTHER = 8
    F_OTHER_MASK = F_OTHER * (FMAX // F_OTHER - 1)

    _Excited = 'E'
    _excited = 'e'
    ISOMER = 'm'
    GROUND = 'g'
    EXCITE = (_Excited, _excited, ISOMER)
    EANY = EMAX - 1
    EXCITED_ANY = '*'
    VOID  = (0, -1, 0, 0, 0)
    VOID_IDX = -1 * FMUL
    VOID_STRING  = '-'
    GAMMA = (0, F_BOSON, 0, 0, +1)
    VAC = (0, 0, 0, 0, 0)
    VAC_STRING = '.'
    VAC_IDX = 0

    E_LEPTON_FLAVOR_MUL = 2
    # lepton number is 2 * (E % 2) - 1, stored a 0, 1 for -1, 1
    # lepton flavors are returned 1-3, stored a 0, 2, 4

    _SPECIAL = {
        0: 'nt',
        1 * EMUL: 'g',
        (F_LEPTON * FMUL) + (+0 * EMUL) + (     +1  * ZMUL) : 'e+',
        (F_LEPTON * FMUL) + (+1 * EMUL) + ((ZMAX-1) * ZMUL) : 'e-',
        (F_LEPTON * FMUL) + (+2 * EMUL) + (     +1  * ZMUL) : 'm+',
        (F_LEPTON * FMUL) + (+3 * EMUL) + ((ZMAX-1) * ZMUL) : 'm-',
        (F_LEPTON * FMUL) + (+4 * EMUL) + (     +1  * ZMUL) : 't+',
        (F_LEPTON * FMUL) + (+5 * EMUL) + ((ZMAX-1) * ZMUL) : 't-',
        (F_LEPTON * FMUL) + (+0 * EMUL) : 'nbe',
        (F_LEPTON * FMUL) + (+1 * EMUL) : 'nue',
        (F_LEPTON * FMUL) + (+2 * EMUL) : 'nbm',
        (F_LEPTON * FMUL) + (+3 * EMUL) : 'num',
        (F_LEPTON * FMUL) + (+4 * EMUL) : 'nbt',
        (F_LEPTON * FMUL) + (+5 * EMUL) : 'nut',
        VOID_IDX : VOID_STRING,
        VAC_IDX : VAC_STRING,
        }
    SPECIAL = {v:k for k,v in _SPECIAL.items()}

    def __init__(self,
                 name = None,
                 Z = None,
                 A = None,
                 E = None,
                 F = None,
                 B = None,
                 idx = None,
                 N = None,
                 element = False,
                 isomer = False,
                 isobar = False,
                 isotone = False,
                 lepton = False,
                 void = True,
                 factory = False):
        """
        Initialize Isotope.

        isomer == True:
            interpret istopes as gs isomers
        """
        if self.__class__ == Ion and factory is False:
            raise NotImplementedError('You should not call Ion constructor directly')

        self.F = self.F_ISOMER if F is None else F
        self.Z = 0 if Z is None else Z
        self.A = 0 if A is None else A
        self.E = 0 if E is None else E
        self.B = 0 if B is None else B
        if N is not None:
            if Z is not None:
                self.A = Z + N
                assert (A == self.A) or (A is None), "Conflict."
            elif A is not None:
                self.Z = A - N
                self.A = A
            else:
                self.F = self.F_ISOTONE
                self.A = N
                assert (self.F == F) or (F is None), "Conflict."
                assert self.E == 0, "Not sure what this would be."
                assert self.A >= 0, "Not sure what this would be."
        if Z is not None:
            if (N is None) and (A is None) and (E is None):
                self.F = self.F_ELEMENT
                assert (self.F == F) or (F is None), "Conflict."
                assert self.E == 0, "Not sure what this would be."
            elif B is None:
                assert self.A != 0,  "Not sure what this would be."
        if A is not None:
            if (Z is None) and (N is None):
                self.F = self.F_ISOBAR
                assert (self.F == F) or (F is None), "Conflict."
                assert self.E == 0, "Not sure what this would be."
                assert self.A > 0, "Not sure what this would be."
        if (E is None) or (E == -1):
            if self.F & self.F_GROUP_MASK == self.F_ISOMER:
                if not isomer:
                    self.F = self.F_ISOTOPE
                self.E == 0
        s = name
        if self._is_ion(s):
            self.B, self.F, self.Z, self.A, self.E = s.tuple()
            # conversion to Ion class
            # or could add keyword
            if not issubclass(type(self), type(s)):
                self.B = 0
                self.F = self.F & self.F_GROUP_MASK
            s = None

        if idx is None:
            try:
                i = int(s)
            except (ValueError, TypeError):
                pass
            else:
                if isobar:
                    self.A = i
                    assert self.Z == 0
                    assert self.E == 0
                    assert self.B == 0
                    self.F = self.F_ISOBAR
                    assert (self.F == F) or (F is None), "Conflict."
                elif isotone:
                    self.A = i
                    assert self.Z == 0
                    assert self.E == 0
                    assert self.B == 0
                    self.F = self.F_ISOTONE
                    assert (self.F == F) or (F is None), "Conflict."
                elif element:
                    self.Z = i
                    assert self.A == 0
                    assert self.E == 0
                    assert self.B == 0
                    self.F = self.F_ELEMENT
                    assert (self.F == F) or (F is None), "Conflict."
                elif lepton:
                    self.Z = i
                    if self.Z == 0 and E == None:
                        self.E = 1
                    else:
                        self.E = E if E is not None else 0
                    assert self.A == 0
                    assert self.B == 0
                    self.F = self.F_LEPTON
                    assert (self.F == F) or (F is None), "Conflict."
                elif 0 < i < len(Elements):
                    self.Z = i
                    self.A = Z if Z is not None else 0
                    self.E = A if A is not None else 0
                    self.F = self.F_ISOMER
                    if A is None:
                        self.F = self.F_ISOTOPE
                    if Z is None:
                        self.F = self.F_ELEMENT
                    assert (self.F == F) or (F is None), "Conflict."
                else:
                    idx = s
                s = None

        if isinstance(s, str):
            try:
                x = tuple(eval(re.sub('[- \.;:/\|]', ',', s)))
                s = tuple(int(i) for i in x)
            except:
                pass
        if isinstance(s, tuple):
            if len(s) == 1:
                self.Z = s[0]
                self.F = self.F_ELEMENT
            if len(s) == 2:
                self.Z = s[0]
                self.A = s[1]
                self.F = self.F_ISOTOPE
            if len(s) == 3:
                self.Z = s[0]
                self.A = s[1]
                self.E = s[2]
                self.F = self.F_ISOMER
            if len(s) == 4:
                self.F = s[0]
                self.Z = s[1]
                self.A = s[2]
                self.E = s[3]
            if len(s) == 5:
                self.B = s[0]
                self.F = s[1]
                self.Z = s[2]
                self.A = s[3]
                self.E = s[4]
            s = None
        if isinstance(s, Ion):
            self.B, self.F, self.Z, self.A, self.E = s.tuple()
        elif s is not None:
            self.B, self.F, self.Z, self.A, self.E = self.ion2bfzae(s, element = element, isomer = isomer)
        if idx is not None:
            self.B, self.F, self.Z, self.A, self.E = self.idx2bfzae(idx)
        self.B = int(self.B)
        self.F = int(self.F)
        self.Z = int(self.Z)
        self.A = int(self.A)
        self.E = int(self.E)
        self.idx = self.bfzae2idx(B = self.B, F = self.F, Z = self.Z, A = self.A, E = self.E)
        # some checks?
        if self.F in (self.F_ISOTOPE, self.F_ISOMER):
            if self.A < self.Z:
                self.idx = self.VOID_IDX
                self.B, self.F, self.Z, self.A, self.E = self.idx2bfzae(self.idx)
        if self.idx == self.VOID_IDX and not void:
            raise AttributeError('Could not create valid isotope from provided input.')
        self.N = self.get_N()
        if isinstance(self, Ion):
            super().__init__()

    @staticmethod
    def _is_ion(ix):
        return  isinstance(ix, Ion)

    @classmethod
    def min_iso(cls, Z):
        if issubclass(type(Z),cls):
            Z = Z.Z
        try:
            Z = int(Z)
        except:
            Z = el2z[Z]
        return cls(Z = Z, N = 0)

    @classmethod
    def max_iso(cls, Z):
        if issubclass(type(Z),cls):
            Z = Z.Z
        try:
            Z = int(Z)
        except:
            Z = el2z[Z]
        return cls(Z = Z, A = cls.AMAX -1)

    @classmethod
    def element_slice(cls, Z):
        return slice(cls.min_iso(Z), cls.max_iso(Z))

    @classmethod
    def bfzae2idx(cls,
                  F = F_ISOTOPE,
                  Z = 0,
                  A = 0,
                  E = 0,
                  B = 0):
        """
        Compute index from F, Z, A, E.
        """
        assert ((F == -1) and (A == 0) and (Z == 0) and (E == 0) and (B == 0)) \
               or ((F >= 0) and (A >= 0) and (E >= 0) and (B >= 0)) \
               or (F >= cls.F_OTHER), "Wrong numbers."
        if Z < 0:
            Z += cls.ZMAX
        return (F * cls.FMUL +
                E * cls.EMUL +
                Z * cls.ZMUL +
                A * cls.AMUL +
                B * cls.BMUL)

    @classmethod
    def idx2bfzae(cls, idx):
        """
        Return things that are defined within basic Ion class.
        """
        aidx = abs(idx)
        f = (aidx // cls.FMUL) % cls.FMAX
        z = (aidx // cls.ZMUL) % cls.ZMAX
        a = (aidx // cls.AMUL) % cls.AMAX
        e = (aidx // cls.EMUL) % cls.EMAX
        b = (aidx // cls.BMUL)

        if z > cls.ZLIM:
            z -= cls.ZMAX
        elif idx == cls.VOID_IDX:
            b, f, z, a, e = cls.VOID
        elif idx == cls.VAC_IDX:
            b, f, z, a, e = cls.VAC
        return b, f, z, a, e

    @classmethod
    def ion2idx(cls,s=''):
        return cls.bfzae2idx(cls.ion2bfzae(s))

    def update_idx(self):
        """
        Update idx assuming B,F,A,Z,E are all current
        """
        self.idx = (self.F * self.FMUL +
                    self.E * self.EMUL +
                    self.Z * self.ZMUL +
                    self.A * self.AMUL +
                    self.B * self.BMUL )

    def get_N(self):
        N = self.A - self.Z
        if self.F in (self.F_BOSON,
                      self.F_ELEMENT,
                      self.F_ISOBAR,
                      self.F_HADRON,
                      self.F_LEPTON):
            N = 0
        return N

    def tuple(self):
        return (self.B, self.F, self.Z, self.A, self.E )

    def dict(self):
        return { 'F': self.F,  'Z': self.Z, 'A': self.A, 'E': self.E, 'B' : self.B }

    def index(self):
        return self.idx

    def index_ZAE(self):
        return (
            self.E * self.EMUL +
            self.Z * self.ZMUL +
            self.A * self.AMUL
            )

    def ZAE(self):
        return (
            self.Z,
            self.A,
            self.E,
            )

    # A,Z, ... should become properties

    def Ye(self):
        if self.Z < 0 or self.A == 0: return -1
        return self.Z / self.A

    def mu(self):
        if self.Z < 0: return -1
        return self.A / (self.Z + 1)

    def mue(self):
        if self.Z <= 0: return -1
        return self.A / self.Z

    def Name(self,
             align = 1,
             width = 0,
             upcase = True):
        """
        Return capitalized ion name; same as name otherwise.
        """
        return self.name(align = align,
                         width = width,
                         upcase = upcase)

    def name(self,
             align = 1,
             width = 0,
             upcase = False):
        """
        Return ion name.

        align [+1]:
            -1: left
             0: center
            +1: right (default)
        width [0]:
            size of field
        upcase [False]:
            capitalize first letter
        """
        s = self._name(upcase = upcase)
        if width > 0:
            if align == -1:
                s = s.ljust(width)
            elif align == 1:
                s = s.rjust(width)
            else:
                s = s.center(width)
        return s

    def element_symbol(self, upcase = True):
        s = ''
        if self.Z >= 0:
            if upcase:
                s = Elements[self.Z]
            else:
                s = elements[self.Z]
        return s

    def LaTeX(self,
              math = False, **kwargs):
        if self.is_element():
            el = Elements[self.Z]
            if math:
                return r"\mathrm{{{:s}}}".format(el)
            else:
                return el
        a = str(self.A)
        el = Elements[self.Z]
        if self.is_isomer():
            kw = {k:kwargs[k] for k in ('m1', 'g') if k in kwargs}
            a += self.isomer_name(self.E, **kw)
        if math:
            return r'^{{{:s}}}\mathrm{{{:s}}}'.format(a, el)
        else:
            return r'$^{{{:s}}}\mathrm{{{:s}}}$'.format(a, el)
        return s

    def NuGrid(self):
        el = Elements[self.Z]
        if self.is_element():
            return el
        A = str(self.A)
        if self.is_isotope():
            return '{}-{}'.format(el, A)
        if self.is_isomer():
            return '{}-{}{}'.format(el, A, self.isomer_name(self.E, m1 = True))
        raise Exception('Cannot convert ion to string.')

    def nugrid5(self):
        el = elements[self.Z]
        if self.is_isotope():
            s = '{:2s}{:3d}'.format(el, self.A)
            if s == "nt  1":
                s = 'NEUT'
            elif s == "h   1":
                s = 'PROT'
            return s
        raise Exception('Cannot convert ion to string.')

    def nugrid5z(self):
        if self.is_isotope():
            return '{:3d} {:5s}'.format(self.Z, self.nugrid5())
        raise Exception('Cannot convert ion to string.')

    def Kepler(self):
        el = elements[self.Z]
        if self.is_element():
            return el
        A = str(self.A)
        if self.is_isotope():
            return '{}{}'.format(el, A)
        if self.is_isomer():
            return '{}{}{}'.format(el, A, self.isomer_name(self.E, m1 = True))
        raise Exception('Cannot convert ion to string.')

    @classmethod
    def isomer_name(cls, E, m1 = False, g = True):
        """
        For isomers use 'g', 'm', 'm2', m3', ...
        """
        if E == 0:
            if g is False:
                return ''
            return cls.GROUND
        elif E == 1 and m1 is False:
            return "{:s}".format(cls.ISOMER)
        elif E == cls.EANY:
            return "{:s}".format(cls.EXCITED_ANY)
        else:
            return "{:s}{:d}".format(cls.ISOMER, E)

    def mpl(self):
        if self.is_element():
            el = Elements[self.Z]
            return r"$\mathsf{{{:}}}$".format(el)
        if self.is_isobar():
            return r"$\mathsf{{{:d}}}$".format(self.A)
        if self.is_isotone():
            return r"$\mathsf{{{:d}}}$".format(self.N)
        if self.is_isotope():
            el = Elements[self.Z]
            a = str(self.A)
            return r"$\mathsf{{^{{{:}\!}}{:}}}$".format(a,el)
        if self.is_isomer():
            el = Elements[self.Z]
            a = str(self.A)
            m = self.isomer_name(self.E)
            return r"$\mathsf{{^{{{:}\!}}{:}^{{{:}\!}}}}$".format(a,el,m)
        return '*'

    def LaTeX_Table(self,
                    math = False,
                    nozero = False,
                    protons = False,
                    manual = False,
                    width = 11,
                    align = 0
                    ):
        a = str(self.A)
        el= Element[self.Z]
        if manual:
            s = '$^{'+a+'}$&'+el
        else:
            s = "\I{"+a+"}{"+el+"}"
        if nozero and idx == VOID_IDX:
            s = '--'
        elif nozero and idx == VAC_IDX:
            s = VAC_STRING
        if align == -1:
            s = s.ljust(width)
        elif align == 1:
            s = s.rjust(width)
        else:
            s = s.center(width)
        if math and not manual: s = '$'+s+'$'
        return s

    def _name(self, upcase = True):
        """
        Convert normal isotopes and isomers to string.

        For isomers use 'g', 'm', 'm2', m3', ...
        """
        s = ''
        if self.Z >= 0:
            if upcase:
                s = Elements[self.Z]
            else:
                s = elements[self.Z]
        if self.F & self.F_GROUP_MASK == self.F_ISOBAR:
            s = 'A:'
        if self.F & self.F_GROUP_MASK == self.F_ISOTONE:
            s = 'N:'
        if self.A != 0 or (self.F & self.F_GROUP_MASK == self.F_ISOTONE):
            s += "{:d}".format(self.A)
        if self.A == 1 and self.Z == 0 and (self.F & self.F_GROUP_MASK == self.F_ISOTOPE):
            s = 'n'
        if self.F & self.F_GROUP_MASK == self.F_ISOMER:
            if self.A == 0 and self.Z == 0:
                if self.E == 1:
                    s = 'g'
                else:
                    s = 'g{:d}'.format(self.E)
            else:
                s += self.isomer_name(self.E)
        if self.F & self.F_GROUP_MASK == self.F_BOSON:
            if self.A == 0 and self.Z == 0:
                if self.E == 1:
                    s = 'g'
                else:
                    s = 'g{:d}'.format(self.E)
            else:
                raise NotImplementedError()
        s = self._SPECIAL.get(self.idx, s)
        return s

    def lepton_number(self):
        """
        return lepton number
        """
        if self.is_lepton:
            return 2 * (self.E % 2) - 1
        return 0

    def lepton_flavor(self):
        """return lepton flavor 1-3, 0 for non-leptons"""
        if self.is_lepton:
            return self.E // 2 + 1
        return 0

    def lepton_flavor_vec(self):
        """return lepton flavour vector [l_e, l_nu, l_tau]"""
        v = [0, 0, 0]
        if self.is_lepton:
            v[self.lepton_flavor() - 1] = self.lepton_number()
        return v

    @classmethod
    def ion2bfzae(cls, s = '', element = False, isomer = False):
        """Interpret string and return component values."""
        s, sz, sa, se = cls._decompose(s, element = element)

        if sz == sa == se == '':
            return cls.VOID

        valid = True
        f = cls.F_ISOMER
        b = 0
        z = None
        try:
            a = int(sa)
        except ValueError:
            a = 0
            if sa != '':
                valid = False
        if sz in ('a=', 'A=', 'A:', 'a:', 'a-', 'A-', 'a', 'A'):
            if a == 0:
               valid = False
            f = cls.F_ISOBAR
            se = 0
            z = 0
        if sz in ('n=', 'N=', 'N:', 'n:', 'n-', 'N-'):
            if a == 0:
               valid = False
            f = cls.F_ISOTONE
            se = 0
            z = 0
        if sz in ('z=', 'Z=', 'Z:', 'Z:', 'z-', 'Z-', 'z', 'Z', 'e', 'E'):
            valid = a > 0
            z = a
            a = 0
            try:
                z = int(sa)
            except ValueError:
                valid = False
            else:
                f = cls.F_ELEMENT
                se = 0
        if z is None:
            try:
                z = Elements.index(sz.strip(',;-=:'))
            except ValueError:
                pass
        if z is None:
            try:
                z = elements.index(sz.lower())
            except ValueError:
                pass
        if z is None:
            try:
                z = z2Name.index(sz.strip('-'))
            except ValueError:
                pass
        if z is None:
            try:
                z = z2name.index(sz.strip('-'))
            except ValueError:
                z = 0
                if sz != '':
                    valid = False
        try:
            e = int(se)
        except ValueError:
            e = 0
            if not isomer:
                f = cls.F_ISOTOPE
            if se != '':
                valid = False
        if sz in ('nt', 'Nt' ) and a == 0:
            f = cls.F_ELEMENT
        elif sz in ('n', 'N') and a > 1:
            z = 7
        elif sz in ('p', 'P') and a > 1:
            z = 15
        elif sz.lower() in ('g', 'gam' ,'gamma') and a > 0:
            f = cls.F_BOSON
            e = a
            a = 0
            valid = True
        elif z == a == 0 and e > 0:
           f = cls.F_BOSON
        else:
            try:
                i = cls.SPECIAL[s]
                b, f, z, a, e = cls.idx2bfzae(i)
                valid = True
            except:
                pass
        if f == cls.F_ISOTOPE and z > 0 and a == 0:
            f = cls.F_ELEMENT
        if f == cls.F_ISOMER and z > 0 and a == 0:
            valid = False
        if not valid:
            b, f, z, a, e = cls.VOID
        return b, f, z, a, e

    def isomer(self, Z = None, A = None, N = None, E = None):
        """
        return self if isomer, isomer if enough info is provided, or VOID otherwise

        """
        if self.is_isomer() and (E is None or isinstance(E, Mapping)):
            if A == N == Z == None:
                return self
            E = self.E
        x = self.isotope(A = A, Z = Z, N = N)
        if x is VOID:
            return ion
        if E == None:
            E = isomermap
        if isinstance(E, Mapping):
            try:
                E = E[x]
            except KeyError:
                # should default behavior be to set E = 0 instead?
                return VOID
        if E is None:
            return VOID
        return Isomer(A = x.A, Z = x.Z, E = E)

    def isotope(self, Z = None, A = None, N = None):
        """
        return isotope if isomer, or self if isotope, or VOID otherwise
        """
        if self.is_isotope() and Z == A == N is None:
            return self

        n = 3 - [Z, A, N].count(None)
        if self.is_isotope() or self.is_isomer():
            if n < 2 and Z is None:
                Z = self.Z
                n += 1
            if n < 2 and A is None:
                A = self.A
                n += 1
        elif self.is_isobar(): # has A
            if n < 2 and A is None:
                A = self.A
        elif self.is_isotone(): # has N
            if n < 2 and N is None:
                N = self.N
        elif self.is_element(): # has Z
            if n < 2 and Z is None:
                Z = self.Z

        # create output
        if Z is None:
            if A is None:
                return VOID
            else:
                if N is None:
                    return VOID
                else:
                    if N > A:
                        return VOID
                    return Isotope(N = N, A = A)
        else:
            if A is None:
                if N is None:
                    return VOID
                else:
                    return Isotope(Z = Z, N = N)
            else:
                if N is None:
                    if Z > A:
                        return VOID
                    return Isotope(Z = Z, A = A)
                else:
                    if N + Z == A:
                        return Isotope(Z = Z, A = A)
        return VOID

    def isobar(self, A = None):
        """
        return isobar if isomer, isotope, or self if isobar, or VOID otherwise
        """
        if A is not None:
            return Isobar(A = A)
        if self.is_isobar():
            return self
        if (self.is_isotope() or self.is_isomer()):
            return Isobar(A = self.A)
        return self.VOID

    def isotone(self, N = None):
        """
        return isotone if isomer, or isotope, or self if isotone, or VOID otherwise
        """
        if N is not None:
            return Isotone(N = N)
        if self.is_isotone():
            return self
        if (self.is_isotope() or self.is_isomer()):
            return Isotone(N = self.N)
        return self.VOID

    def element(self, Z = None):
        """
        return element if isomer, or isotope, or self if element, or VOID otherwise
        """
        if Z is not None:
            return Element(Z = Z)
        if self.is_element():
            return self
        if (self.is_isotope() or self.is_isomer()):
            return Element(Z = self.Z)
        return self.VOID

    _translate_nugrid = {
        'cd*15' : 'cd115m',
        'lu*76' : 'lu115m',
        'tag80' : 'ta180g',
        'prot'  : 'pn1',
        'neut'  : 'nt1',
        'OOOOO' : 'g',
        }

    _translate_reaclib = {
        'al-6'  : 'al26g',
        'al*6'  : 'al26m',
        }

    _translate_extra = {
        'alpha' : 'he4',
        'gamma' : 'g',
        }

    _translate = {}
    _translate.update(_translate_nugrid)
    _translate.update(_translate_reaclib)
    _translate.update(_translate_extra)

    _html = re.compile('<sup>([0-9mg*]+)</sup>([a-zA-Z]+)')

    @classmethod
    def _decompose(cls,
                   s = '',
                   element = False):
        """
        This routine TRIES to decompose string in
        1) element symbol
        2) mass number
        3) excited state number
        All are returned as strings.

        The optional parameter 'element' enforces that
        'p' is taken as *element* potassium (not H1)
        'n' is taken as *element* nitrogen (not neutron)
        """

        s = s.strip()

        x = cls._html.findall(s)
        if len(x) > 0:
            s = ''.join(x[0][::-1])

        s = cls._translate.get(s.lower(), s)

        name = s.strip()
        n = len(name)
        el = ''
        a = ''
        e = ''

        # get numbers
        n = re.findall("\d+", name)

        # get strings
        cx = re.findall("\D+", name)

        c = []
        for x in cx:
            xx = x.split('-')
            cy = [y for y in xx if y !=  '']
            c += cy
        if len(c) == 2:
            if c[0] in ('m', 'g'):
                c = c[::-1]
            if c[0][0] == '*':
                c = c[::-1]
        if len(n) > 0: a = n[0]
        if len(n) > 1: e = n[1]
        if len(n) > 2: raise ValueError("Can't understand isotope '{}'.".format(s))
        if len(c) > 0: el = c[0]
        if len(el) > 0:
            if el[-1] in cls.EXCITE and len(c) == 1 and len(n) == 2:
                c.append(el[-1])
                el = el[:-1]
        if len(c) == 2 and c == ['(', ')']:
            if len(n) == 1:
                a = n[0]
                el = 'Z='
                e = ''
                c = []
                n = []
            else:
                return (s,) + ('',)*3
        if len(c) == 2:
            if c[1] in ('g', 'G'):
                e = '0'
                if len(n) > 1:
                    return (s,) + ('',)*3
            elif c[1] in ('m', 'M') and len(n) == 1:
                e = '1'
            elif c[1][0] == '*' and len(n) == 1:
                e = str(len(c[1]))
                assert c[1].count('*') == len(c[1])
                if e == '1':
                    e = str(cls.EANY)
            if not c[1] in ('m', 'g', 'M', 'G') and not c[1][0] == '*':
                return (s,) + ('',)*3

        if len(c) == 1 and c[0][-1] == '*':
            e = 0
            while c[0][-1] == '*':
                c[0] = c[0][:-1]
                e += 1
            assert e == 1
            e = str(e)
            el = c[0]

        if len(c) == 1 and c[0][0] == '*':
            e = 0
            while c[0][0] == '*':
                c[0] = c[0][1:]
                e += 1
            assert e == 1
            e = str(e)
            el = c[0]

        if s == 'a' and a == '':
            el = 'He'
            a = '4'
        # this is a possible conflict with potassium
        elif (element) and s == 'p':
            el = 'P'
        elif s == 'p':
            el = 'H'
            a = '1'
        elif el in ('p', 'pn') and a == '1':
            el = 'H'
        elif s == 'pn':
            el = 'H'
            a = ''
        elif el in ('d', 'D'):
            el = 'H'
            if not a in ('', '2'):
                raise AttributeError('"d" already implies mass; if supplied needs to be "2".')
            a = '2'
        elif el in ('t','T'):
            el = 'H'
            if not a in ('', '3'):
                raise AttributeError('"t" already implies mass; if supplied needs to be "3"')
            a = '3'
        elif (element) and s == 'n':
            el = 'N'
        elif s == 'n':
            el = 'nt'
            a = '1'
        elif el in ('n', 'nt') and a == '1':
            el = 'nt'
        elif s in ('g', 'G'):
            el = ''
            a = ''
            e = '1'
        elif (s.lower() in ('e-', 'b-', 'bd', 'pc')):
            s = el = 'e-'
        elif ((s.lower() in ('e+', 'b+', 'ec'))
            or ((not element) and (s.lower() == 'pd'))):
            s = el = 'e+'
        elif ((not element) and (s.lower() == 'ps')):
            s = 'h1'
            a = '1'
            el = 'h'
        elif ((not element) and (s.lower() == 'ns')):
            s = 'nt1'
            a = '1'
            el = 'nt'
        el = el.strip()
#        if len(el) == 2 and el(2)
        a = a.strip()
        e = e.strip()
        return s, el, a, e

    def deexcite(self):
        self.E = 0
        self.idx = self.bfzae2idx(F = self.F,
                                  Z = self.Z,
                                  A = self.A,
                                  E = self.E,
                                  B = self.B)

    def deexcited(self):
        return ion(
            F = self.F,
            A = self.A,
            Z = self.Z,
            B = self.B,
            )

    def __str__(self):
        return self._name(upcase = True)

    def __format__(self, format_spec):
        return '{:{}}'.format(self.__str__(), format_spec)

    def nugridpy(self):
        # need to update for isotopes
        return '{}-{:d}'.format(
            Elements[self.Z],
            self.A)

    _custom_add = False

    def _add(self, x, sign1 = +1, sign2 = +1):
        if np.shape(x) is not ():
            return NotImplemented
        if not self._is_ion(x):
            y = ion(x)
            if y.idx == self.VOID_IDX:
                raise AttributeError('Cannot convert {!r} to {!r}.'.format(x, Ion.__class__.__name__))
            else:
                x = y
        if x._custom_add:
            return NotImplemented
        A = sign1 * self.A + sign2 * x.A
        Z = sign1 * self.Z + sign2 * x.Z
        E = sign1 * self.E + sign2 * x.E
        if not (((self.F & self.F_GROUP_MASK == self.F_ISOMER) and x.is_photon()) or
                (x.F    & self.F_GROUP_MASK == self.F_ISOMER) and x.is_photon()):
            E = None
        return ion(
            Z = Z,
            A = A,
            E = E,
            )

    def __add__(self, x):
        return self._add(x)

    __radd__ = __add__

    def __sub__(self, x):
        return self._add(x, +1, -1)

    def __rsub__(self, x):
        return self._add(x, -1, +1)

    def __mul__(self, x):
        if np.shape(x) is not ():
             return NotImplemented
        return self.__class__(
            Z = self.Z * x,
            A = self.A * x,
            E = self.E,
            F = self.F,
            )
    __rmul__ = __mul__

    def __neg__(self):
        return self.__class__(
            Z = - self.Z,
            A = - self.A,
            E =   self.E,
            F =   self.F,
            B =   self.B,
            )

    def __bool__(self):
        return self.idx != self.VOID_IDX

    def __hash__(self):
        return self.idx

    def __lt__(self, x):
        if np.shape(x) is not ():
            return NotImplemented
        if not isinstance(x, Ion):
            x = self.factory(x)
        return self.idx < x.idx

    def __eq__(self, x):
        if np.shape(x) is not ():
            return NotImplemented
        if not isinstance(x, Ion):
            x = self.factory(x)
        return self.idx == x.idx

    def __len__(self):
        return 1

    def __repr__(self):
        base = "('{:s}')".format(self.Name())
        #     if self.is_isomer():
        #         return 'Isomer' + base
        #     if self.is_isotope():
        #         return 'Isotope' + base
        #     if self.is_isobar():
        #         return 'Isobar' + base
        #     if self.is_isotine():
        #         return 'Isotone' + base
        #     if self.is_element():
        #         return 'Element' + base
        #     if self.is_lepton():
        #         return 'Lepton' + base
        #     if self.is_hadron():
        #         return 'Hadron' + base
        return self.__class__.__name__ + base

    def __getstate__(self):
        # need to add version number
        # should save index instead
#        print('xxx_save')
        return self.tuple()

    def __setstate__(self,x):
        # need to add check of version number
#        print('xxx_load')
        self.__init__(x)
# this does not appear to be called for a new-style class ever.
    # def __getinitargs__(self):
    #     print 'xxx_old'
    #     return (self.tuple(),)
# the following could work as well, but does appear less flexible
    # def __getnewargs__(self):
    #     print 'xxx_new'
    #     return (self.tuple(),)
    def element_name(self):
        return z2name[self.Z]

    def is_nucleus(self):
        return self.is_isomer() or self.is_isotope()
    def is_lepton(self):
        return self.F & self.F_GROUP_MASK == self.F_LEPTON and self.F != -1
    def is_isomer(self):
        return self.F & self.F_GROUP_MASK == self.F_ISOMER
    def is_isotope(self):
        return self.F & self.F_GROUP_MASK == self.F_ISOTOPE
    def is_isobar(self):
        return self.F & self.F_GROUP_MASK == self.F_ISOBAR
    def is_isotone(self):
        return self.F & self.F_GROUP_MASK == self.F_ISOTONE
    def is_element(self):
        return self.F & self.F_GROUP_MASK == self.F_ELEMENT
    def is_hadron(self):
        return self.F & self.F_GROUP_MASK == self.F_HADRON
    def is_boson(self):
        return self.F & self.F_GROUP_MASK == self.F_BOSON
    def is_photon(self):
        return ((self.F & self.F_GROUP_MASK == self.F_BOSON)
                and (self.Z == 0)
                and (self.A == 0)
                and (self.E > 0)
                )
    def is_void(self):
        return self.F == -1
    def is_vac(self):
        return self.F == 0

    def type(self):
        if self.is_lepton():
            return 'lepton'
        if self.is_isotope():
            return 'isotope'
        if self.is_isobar():
            return 'isobar'
        if self.is_isotone():
            return 'isotone'
        if self.is_element():
            return 'element'
        if self.is_isomer():
            return 'isomer'
        if self.is_hadron():
            return 'hadron'
        if self.is_nucleon():
            return 'nucleon'
        if self.is_photon():
            return 'photon'
        if self.is_boson():
            return 'boson'
        if self.is_void():
            return 'void'
        if self.is_vac():
            return 'vac'
        return 'unkown'

    def anti_particle(self):
        return self.__class__(
            A = self.A,
            Z = -self.Z,
            F = self.F,
            B = self.B,
            )

    def R(self,
          I = None,
          O = None):
        """
        Reaction rate function.
        If only one parameter is given it is interpreted as 'outgoing'
        If no parameter is provided, just a deep copy of self is returned
        """
        if O is None:
            if I is None:
                return self()
            return self(None, I)
        return self(I, O)

    def __call__(self, *args):
        """
        implement reaction semantics

        if only one parameter is supplied, assume it is the output
        channel

        if more than one is supplied, assume the first is the in
        channel, the other are the out channel

        in/out channels may be iterables

        Examples
        ion('c12')('p','g') --> Isotope('N13')

        ion('c12')('p') -->

        """
        result = self
        if len(args) == 1:
            if np.isscalar(args[0]) or args[0] is None:
                result -= args[0]
            else:
                for i in args[0]:
                    result -= i
            return result
        if np.isscalar(args[0]) or args[0] is None:
            result += args[0]
        else:
            for i in args[0]:
                result += i
        for i in args[1:]:
            if np.isscalar(i) or i is None:
                result -= i
            else:
                for j in i:
                    result -= j
        return result

ufunc_A = lambda y: np.array(np.frompyfunc(lambda x: x.A, 1, 1)(y), dtype = np.int)
ufunc_Z = lambda y: np.array(np.frompyfunc(lambda x: x.Z, 1, 1)(y), dtype = np.int)
ufunc_N = lambda y: np.array(np.frompyfunc(lambda x: x.N, 1, 1)(y), dtype = np.int)
ufunc_E = lambda y: np.array(np.frompyfunc(lambda x: x.E, 1, 1)(y), dtype = np.int)
ufunc_F = lambda y: np.array(np.frompyfunc(lambda x: x.F, 1, 1)(y), dtype = np.int)
ufunc_B = np.frompyfunc(lambda x: x.B, 1, 1) # these may be larger than int64

ufunc_isomer  = lambda y: np.array(np.frompyfunc(lambda x: x.isomer() , 1, 1)(y), dtype = np.object)
ufunc_isotope = lambda y: np.array(np.frompyfunc(lambda x: x.isotope(), 1, 1)(y), dtype = np.object)
ufunc_element = lambda y: np.array(np.frompyfunc(lambda x: x.element(), 1, 1)(y), dtype = np.object)
ufunc_isobar  = lambda y: np.array(np.frompyfunc(lambda x: x.isobar() , 1, 1)(y), dtype = np.object)
ufunc_isotone = lambda y: np.array(np.frompyfunc(lambda x: x.isotone(), 1, 1)(y), dtype = np.object)

ufunc_element_name   = lambda y: np.array(np.frompyfunc(lambda x: x.element_name()               , 1, 1)(y))
ufunc_element_symbol = lambda y: np.array(np.frompyfunc(lambda x: x.element_symbol(upase = True) , 1, 1)(y))
ufunc_element_sym_lc = lambda y: np.array(np.frompyfunc(lambda x: x.element_symbol(upase = False), 1, 1)(y))

ufunc_isomer_idx  = lambda y: np.array(np.frompyfunc(lambda x: x.isomer().idx , 1, 1)(y), dtype = np.int)
ufunc_isotope_idx = lambda y: np.array(np.frompyfunc(lambda x: x.isotope().idx, 1, 1)(y), dtype = np.int)
ufunc_element_idx = lambda y: np.array(np.frompyfunc(lambda x: x.element().idx, 1, 1)(y), dtype = np.int)
ufunc_isobar_idx  = lambda y: np.array(np.frompyfunc(lambda x: x.isobar().idx , 1, 1)(y), dtype = np.int)
ufunc_isotone_idx = lambda y: np.array(np.frompyfunc(lambda x: x.isotone().idx, 1, 1)(y), dtype = np.int)

ufunc_is_ion  = lambda y: np.array(np.frompyfunc(lambda x: isinstance(x, Ion), 1, 1)(y), dtype = np.bool)

ufunc_is_lepton  = lambda y: np.array(np.frompyfunc(lambda x: x.is_lepton() , 1, 1)(y), dtype = np.bool)
ufunc_is_isotope = lambda y: np.array(np.frompyfunc(lambda x: x.is_isotope(), 1, 1)(y), dtype = np.bool)
ufunc_is_isobar  = lambda y: np.array(np.frompyfunc(lambda x: x.is_isobar() , 1, 1)(y), dtype = np.bool)
ufunc_is_isotone = lambda y: np.array(np.frompyfunc(lambda x: x.is_isotone(), 1, 1)(y), dtype = np.bool)
ufunc_is_element = lambda y: np.array(np.frompyfunc(lambda x: x.is_element(), 1, 1)(y), dtype = np.bool)
ufunc_is_isomer  = lambda y: np.array(np.frompyfunc(lambda x: x.is_isomer() , 1, 1)(y), dtype = np.bool)
ufunc_is_hadron  = lambda y: np.array(np.frompyfunc(lambda x: x.is_hadron() , 1, 1)(y), dtype = np.bool)
ufunc_is_nucleus = lambda y: np.array(np.frompyfunc(lambda x: x.is_nucleus(), 1, 1)(y), dtype = np.bool)
ufunc_is_photon  = lambda y: np.array(np.frompyfunc(lambda x: x.is_photon() , 1, 1)(y), dtype = np.bool)
ufunc_is_void    = lambda y: np.array(np.frompyfunc(lambda x: x.is_void()   , 1, 1)(y), dtype = np.bool)
ufunc_is_vac     = lambda y: np.array(np.frompyfunc(lambda x: x.is_vac()    , 1, 1)(y), dtype = np.bool)

ufunc_type = lambda y: np.array(np.frompyfunc(lambda x: x.type, 1, 1)(y), dtype = np.str)

ufunc_idx  = lambda y: np.array(np.frompyfunc(lambda x: x.idx, 1, 1)(y), dtype = np.int)

ufunc_ion  = lambda y: np.array(np.frompyfunc(lambda x: ion(x), 1, 1)(y), dtype = np.object)
ufunc_ion_from_idx  = lambda y: np.array(np.frompyfunc(lambda x: ion(idx = x), 1, 1)(y), dtype = np.object)

# some convenience functions
def ionarr(arr):
    if isinstance(arr, str):
        arr = re.split('[^-0-9a-zA-Z*]+', arr)
    if not isinstance(arr, Iterable):
        arr = (arr,)
    return np.array([ion(a) for a in arr], dtype = np.object)

# the following classes need to be fleshed out and used consequnetly
# in particular, maybe some filter of constructor can be done?
class Element(Ion):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert self.is_element()

class Isotone(Ion):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert self.is_isotone()

class Isobar(Ion):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert self.is_isobar()

class Isotope(Ion):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert self.is_isotope()

class Isomer(Ion):
    def __init__(self, *args, **kwargs):
        kwargs = kwargs.copy()
        kwargs['isomer'] = True
        super().__init__(*args, **kwargs)
        assert self.is_isomer()

class Photon(Ion):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert self.is_photon()

class Lepton(Ion):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert self.is_lepton()

class Hadron(Ion):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert self.is_hadron()

class Void(Ion):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert self.is_void()

class Vac(Ion):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert self.is_vac()

# creation functions/factories
def element(Z):
    return Element(Z = Z)
def isotone(N):
    return Isotone(N = N)
def isobar(A):
    return Isobar(A = A)
def isotope(A, Z):
    return Isotope(A = A, Z = Z)
def isomer(A, Z, E = 0):
    return Isomer(A = A, Z = Z, E = E)
def photon(E = 1):
    return Photon(E = E)
def lepton(Z = 0, E = 1):
    return Lepton(Z = Z, E = E, lepton = True)
def hadron(Z = 0, E = 1):
    return Hadron(Z = Z, E = E, F = Ion.F_HARDON)
def void():
    return Void(Ion.VOID)

class IonFactory():
    """
    Factory function class.

    This one should be used to allow for later refactoring in derived
    classes for different ion types.
    """
    def __call__(self, *args, **kwargs):
        if len(args) == 1 and len(kwargs) == 0 and isinstance(args[0], Ion):
            return args[0]
        ix = Ion(*args, factory = True, **kwargs)
        if ix.is_isomer():
            return Isomer(ix)
        if ix.is_isotope():
            return Isotope(ix)
        if ix.is_element():
            return Element(ix)
        if ix.is_isobar():
            return Isobar(ix)
        if ix.is_isotone():
            return Isotone(ix)
        if ix.is_photon():
            return Photon(ix)
        if ix.is_lepton():
            return Lepton(ix)
        if ix.is_hadron():
            return Hadron(ix)
        if ix.is_void():
            return Void(ix)
        if ix.is_vac():
            return Vac(ix)
        debug = kwargs.get('debug', True)
        if debug:
            raise Exception('Unkown Ion type.')
        return ix
    def __getattr__(self, attr):
        if not isinstance(attr, Ion):
            x = ion(attr)
        else:
            x = attr
        if x is not Ion.VOID:
            return self(x)
        raise AttributeError(attr)
ion = IonFactory()
Ion.factory = ion

VAC = ion(Ion.VAC)
GAMMA = ion(Ion.GAMMA)
VOID  = ion(Ion.VOID)
NEUTRON = ion(Z=0,A=1)
PROTON = ion(Z=1,A=1)
ALPHA = ion(Z=2,A=4)

class IsomerMap(collections.defaultdict):
    """
    define map that can be used as 'E' paramater in ion.isomer function

    TODO - add more isotopes, provide function to load from file
    """
    def __init__(self, default = lambda: 0, isomap = None):
        super().__init__(default)
        if map is not None:
            self.update(isomap)
    def __getitem__(self, item):
        if not isinstance(item, Ion):
            item = ion(item)
        return super().__getitem__(item)
    def __setitem__(self, item, value):
        if not isinstance(item, Ion):
            item = ion(item)
        super().__setitem__(item, value)
isomermap = IsomerMap(isomap = {Isotope('ta180') : 1})

# registry of other bits [set, but not yet used except to check for duplicates]
other_bits_register = {}
def register_other_bits(bit, name = None, class_ = None, overwrite = True):
    if isinstance(bit, type):
        class_ = bit
        name = class_.__name__
        bit = class_.__dict__['F_OTHER_' + name.upper()]

    if name is not None and class_ is None:
        name, class_ = class_, name
    if name is None:
        name = class_.__name__

    # this may need to become more sophisticated to allow re-compilation
    if not overwrite:
        assert bit not in other_bits_register, 'bit already exisits'

    other_bits_register[bit] = (name, class_)
