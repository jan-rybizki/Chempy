"""
Module for reaction based on Ion class

TODO,

"""
import copy
import re
import io
import os
import time

from collections.abc import Iterable

import numpy as np
from numpy.linalg import matrix_power
from scipy.linalg import solve_banded

import isotope
from logged import Logged
from utils import stuple, project, CachedAttribute
from collections import Iterable
from version2human import version2human
from utils import cachedmethod

import physconst

class DuplicateIons(Exception):
    """
    Exception raise when ions list or sets illegally contian duplicate ions
    """
    def __init__(self, ions = None):
        if ions is None:
            super().__init__('Duplicate ions.')
        else:
            super().__init__('Duplicate ions: {}'.format(ions))

class IonList(object):
    """
    Unsorted list of ions - fixed positions.

    TODO - Maybe not allow KepIon and Ion mixed - all things with extra flags.

    TODO - add some of AbuSet functionallity, e.g., return all
    Elements, maybe maping matrices (move functionallity here?)
    """

    def __init__(self, ions = None, duplicates = None):
        """
        duplicates:
          True
          False - drop
          None - raise exception
        """
        self._ions = np.array((), dtype = np.object)
        self.duplicates = duplicates
        self.add(ions)

    def add(self, ions = None):
        """
        Add ion(s) to list.
        """
        if ions is None:
            return
        if isinstance(ions, self.__class__):
            self._ions = np.append(self._ions, ions)
            return
        ions = np.atleast_1d(ions)
        ions = np.array([i if isinstance(i, isotope.Ion) else isotope.ion(i) for i in ions])
        if self.duplicates is True:
            self._ions = np.append(self._ions, ions)
        elif self.duplicates is False:
            _, ii = np.unique(ions, return_index = True)
            ions = ions[sorted(ii)]
            ii = np.in1d(ions, self._ions)
            self._ions = np.append(self._ions, ions[~ii])
        elif self.duplicates is None:
            ions = np.append(self._ions, ions)
            ix, nx = np.unique(ions, return_counts = True)
            ix = ix[nx > 1]
            if len(ix) > 0:
                raise DuplicateIons(ix)
            self._ions = ions
        else:
            raise Exception('Invalid mode for "duplicates": {}'.format(self.duplicates))

    def index(self, ion):
        if not isinstance(ion, isotope.Ion):
            ion = isotope.ion(ion)
        i = np.argwhere(self._ions == ion)
        if len(i) == 0:
            return -1
        else:
            return i.flatten()[0]
    def __eq__(self, other):
        # assert isinstance(other, IonList)
        if len(self) != len(other):
            return False
        for i,j in zip(self, other):
            if not i == j:
                return False
        return True
    def __ne__(self, other):
        return not self.__eq__(other)

    def ions(self):
        return copy.deepcopy(self._ions)

    def copy(self):
        new = self.__class__()
        new._ions = self._ions.copy()
        return new

    # copy from IonList:
    def __str__(self):
        return "[{}]".format(','.join((str(i) for i in self._ions)))
    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,str(self))
    def __getitem__(self, index):
        return self._ions[index]
    def __len__(self):
        """
        Return number of isotopes.
        """
        return len(self._ions)
    def __iter__(self):
        for ion in self._ions:
            yield ion
    def __add__(self, other):
        new = self.copy()
        new.add(other)
        return new

    @property
    def A(self):
        return isotope.ufunc_A(self._ions)

    @property
    def Z(self):
        return isotope.ufunc_Z(self._ions)

    @property
    def N(self):
        return isotope.ufunc_N(self._ions)

class IonSet(object):
    """
    Provide sorted list of isotopes.

    TODO: efficient additions of arrays of ions.

    I suppose in principle we do want :set: properties but it also is
    supposed to be sorted by index at all times if _sort is set.

    Maybe there need be 2 kinds to allow interface with ionmap:
    ionmap also allows to have elements, etc., mixed.

    TODO - Maybe not KepIon and Ion mixed - things with extra flags.

    TODO - Do not allow duplicates

    NOTE - I think there should be only one, IonSet or IonList
    """
    def __init__(self, *args, **kwargs):
        self._ions = np.array((), dtype = np.object)
        self._sort = kwargs.get('sort', True)
        self._type = kwargs.get('type', None)
        assert self._type is None or not issubclass(self._type, isotope.Ion), 'type needs to be subclass of Ion'
        self.add(args)

    def copy(self):
        """
        return a copy
        """
        return copy.copy(self)

    def add(self, ions = None):
        """
        Add ion(s) to list.

        Ions need to be same "type"
        """
        if isinstance(ions, str):
            ions = tuple(isotope.ion(i) for i in re.split('[\'",;\s]+', ions))
        if len(ions) == 1:
            try:
                ions = ions[0]
            except:
                pass
        if (ions is not None):
            if np.isscalar(ions):
                self._add_one(ions)
            elif len(ions) == 1:
                while isinstance(ions, Iterable) and not isinstance(ions, str):
                    ions = ions[0]
                self._add_one(ions)
            else:
                for ion in ions:
                    self._add_one(ion, update = False)
                self._update()

    def _add_one(self, ix, update = True):
        """
        add one ion and check compatibility
        """
        if self._type is None:
            self._type = type(ix)
            if not issubclass(self._type, isotope.Ion):
                self._type = type(isotope.ion(ix))
        if isinstance(ix, isotope.Ion) and not isinstance(ix, self._type):
            raise TypeError("All ions need compatible type: {:s}, {:s}".format(
                    str(self._type),str(type(ix))))
        if not isinstance(ix, self._type):
            ix = self._type(ix)
        self._ions = np.append(self._ions, ix)
        if update:
            self._update()

    def _update(self):
        if self._sort:
            self._ions.sort()
        self._ions = np.unique(self._ions)

    def __str__(self):
        return "[{}]".format(','.join((str(i) for i in self._ions)))
    def __repr__(self):
        return "{}({})".format(self.__class__.__name__,str(self))
    def __getitem__(self, index):
        return self._ions[index]
    def __len__(self):
        """
        Return number of isotopes.
        """
        return len(self._ions)
    def __iter__(self):
        i=0
        while i < self.__len__():
            yield self._ions[i]
            i+=1
    def __add__(self, other):
        new = self.copy()
        new.add(other)
        return new

class Abu(object):
    """
    Ion plus single abundance.  Provide some general functionallity.

    TODO - add most of AbuSet functionallity

    TODO  Contain mass information?
    """
    def __init__(self, ion, abu = 0):
        self.ion = ion
        self.abu = np.float64(abu)
    def __str__(self):
        return "{!s}:{:>12.5f}".format(self.ion,self.abu)
    def X(self):
        return self.abu
    def Y(self):
        return self.abu / self.ion.A

class AbuDist(object):
    """
    Ion and distribution information, e.g., structure.

    TDOD - implement
    TDOD - add Abu(Set) and AbuData functionallities
    """
    # should be initialized with Ion and abundace vector,
    # mass and radius would be good.
    pass


class AbuData(object):
    """
    n-D abu data plus extra info.
    ions are in shape[0]

    Probably base class (if not identical) to starfit data.
    """
    def __init__(self,
                 data = None,
                 ions = None,
                 molfrac = None,
                 #fields = None,
                 #formats = None,
                 #values = None,
                 #info = None,
                 ):

        self.ions = IonList(ions)
        self.data = data
        self.molfrac = molfrac

        #self.fields = fields
        #self.formats = formats
        #self.values = values

    def __len__(self):
        return np.product(self.data.shape[:-1])

    def abu(self, molfrac = None):
        """
        abundace
        """
        if molfrac is None:
            molfrac = self.molfrac

        if molfrac == self.molfrac:
            yfac = np.ones(self.data.shape[1])
        else:
            yfac = isotope.ufunc_A(self.ions)
        if molfrac and not self.molfrac:
            yfac = 1 / yfac

        value = np.ndarray(self.data.shape, dtype = np.float64)
        value[...,:] = self.data[...,:] * yfac
        return value

    def updated(self,
                data = None,
                ions = None,
                molfrac = None,
                ):
        """
        return updated copy of self
        """
        new = copy.deepcopy(self)
        new.data    = data.copy()
        new.ions    = IonList(ions.copy())
        if molfrac is not None:
            new.molfrac = molfrac
        return new

    def decayed(self, decay = None, **decpar):
        """
        Return copy of self with yields replaced by decayed yields.
        """
        from ionmap import Decay

        return Decay.Map(self, **decpar)

    def cleaned(self, threshold = 1.e-25, debug = False):
        maxima = np.max(self.data, axis = 0)
        ii, = np.where(maxima > threshold)
        save_data = self.data
        save_ions = self.ions
        self.data = self.data[:, ii]
        self.ions = self.ions.__class__(self.ions[ii])
        if debug:
            print('[{}] removing {}.'.format(
                self.__class__.__name__,
                ', '.join(str(i) for i in (set(save_ions) - set(self.ions))
                          )
                ))
        new = copy.deepcopy(self)
        self.data = save_data
        self.ions = save_ions
        return new

    @property
    def A(self):
        return self.ions.A

    @property
    def Z(self):
        return self.ions.Z

    @property
    def N(self):
        return self.ions.N

class AbuDump(AbuData):
    """
    KEPLER BURN abundance set.  To hold KEPLER ppnb data (plus wind)
    """
    def __init__(self,
                 data = None,
                 ions = None,
                 windb = None,
                 molfrac = None,
                 bottom = None,
                 xm = None,
                 zm = None,
                 ):
        """
        xm is mass of zone in g
        ym is mass coordinate of out zone boundary

        if there is n star data, it should be in zones 1..n zone 0 is
        the inner boundady condition (could become having non-zero data
        in accretion probelms); zone n+1 holds the wind data

        abundace data should be mol fractions or mass fractions

        data on import has format [isotopes, zones+2], but is stored
        internally as [zones+2, isotopes]

        molfrac determines whether the input is in mol fractions;
        wise to set if known; otherwise the code tries automatic determination,
        which may fail.
        """
        if isinstance(data, AbuData):
            self.data = data.data.copy()
            self.ions = data.ions.copy()
            self.molfrac = data.molfrac
            assert ions == None
            assert molfrac == None
            if bottom is None:
                bottom = 1
            self.i0 = bottom
            self.i1 = data.data.shape[1]
            #check for wind - if abu in wind is 0, reduce i1 by 1
            if np.sum(self.data[-1,:]) == 0:
                self.has_wind = False
                self.i1 -= 1
            else:
                self.has_wind = True
        else:
            self.data = data.copy().transpose()
            self.ions = IonList(ions)
            if molfrac is None:
                molfrac = (1 - np.sum(self.data[1,:])) > 0.1
            self.molfrac = molfrac
            self.i0 = 1
            self.i1 = self.data.shape[0] - 1
            if windb is not None:
                self.data[-1,:] = windb
                self.i1 += 1
                self.has_wind = True
            else:
                self.has_wind = False
        self.xm = xm
        self.zm = zm

    def __iter__(self, **kw):
        for ion in self.ions:
            yield (ion, self.ion_abu(ion, **kw))

    def ion_abu(self,
                ion,
                molfrac = False,
                missing = np.nan):
        """
        Return isotope abundace
        """
        if not isinstance(ion, isotope.Ion):
            ion = isotope.ion(ion)
            if ion == isotope.VOID:
                value = np.ndarray(self.data.shape[0], dtype = np.float64)
                value.fill(missing)
                return value
            # return self.ion_abu(ion, molfrac = molfrac, missing = missing)
        yfac = max(1, ion.A)
        if molfrac == self.molfrac:
            yfac = 1
        if molfrac and not self.molfrac:
            yfac = 1 / yfac
        value = np.ndarray(self.data.shape[0], dtype = np.float64)
        i = self.ions.index(ion)
        assert i != -1, 'ion {} not in data set'.format(ion)
        value[self.i0:self.i1] = self.data[self.i0:self.i1,i] * yfac

        value[:self.i0] = missing
        value[self.i1:] = missing
        return value

    @cachedmethod
    def __call__(self, *args, **kw):
        return self.ion_abu(*args, **kw)

    @cachedmethod
    def __getitem__(self, ion):
        return self.ion_abu(ion)

    def __getattr__(self, attr):
        if attr in self.ions:
            return self.ion_abu(attr)
        return super().__getattribute__(attr)

    def abu(self,
            molfrac = None,
            missing = np.nan,
            wind = True,
            bottom = True):
        """
        abundace
        """
        i1 = self.data.shape[0] - ( 0 if wind else 1)
        i0 = 0 if bottom else 1

        if molfrac is None:
            molfrac = self.molfrac

        if molfrac == self.molfrac:
            yfac = np.ones(self.data.shape[1])
        else:
            yfac = isotope.ufunc_A(self.ions)
        if molfrac and not self.molfrac:
            yfac = 1 / yfac
        value = np.ndarray(self.data.shape, dtype = np.float64)
        value[self.i0:self.i1,:] = self.data[self.i0:self.i1,:] * yfac[np.newaxis,:]
        value[:self.i0,:] = missing
        value[self.i1:,:] = missing
        return value[i0:i1]

    def project(self,
                output = 'massfrac',
                xm = None,
                zones = None,
                ):
        """
        Return projected yields as AbusSet object.

        Output = massfrac | molfrac | g | mol | Msun
        """
        if xm is None:
            try:
                xm = self.xm
            except AttributeError:
                raise Exception('Need to have zone mass defined.')
        assert xm.shape[0] == self.data.shape[0]

        if zones is None:
            # exclude "phony" zones but leave in wind
            zones = slice(self.i0, self.i1)

        norm = 1
        if output == 'massfrac':
            molfrac = False
            norm = 1 / np.sum(xm[zones])
        elif output == 'molfrac':
            molfrac = True
            norm = 1 / np.sum(xm[zones])
        elif output == 'mol':
            molfrac = True
        elif output == 'g':
            molfrac = False
        elif output == 'Msun':
            molfrac = False
            norm = 1 / physconst.XMSUN

        if molfrac is None:
            molfrac = self.molfrac

        if molfrac == self.molfrac:
            yfac = np.ones(self.data.shape[1])
        else:
            yfac = isotope.ufunc_A(self.ions)
        if molfrac and not self.molfrac:
            yfac = 1 / yfac

        value = np.tensordot(self.data[zones, :], xm[zones], axes=(0, 0)) * yfac * norm
        return AbuSet(
            iso = self.ions,
            abu = value,
            normalize = False,
            unit = '?')

    def mix(self,
            mbox,
            iterations = 4,
            xm = None,
            zones = None,
            ):
        """
        Mix abundance dump and return new mixed object.
        """
        if xm is None:
            try:
                xm = self.xm
            except AttributeError:
                raise Exception('Need to have zone mass defined.')
        assert xm.shape[0] == self.data.shape[0], 'require same number of zones in xm and data set'

        if zones is None:
            # exclude "phony" and surface zones
            if self.has_wind:
                zones = slice(self.i0, self.i1 - 1)
            else:
                zones = slice(self.i0, self.i1)

        # convert to list of zones
        if isinstance(zones, slice):
            zones = np.array(range(*zones.indices(len(xm))))
        if not isinstance(zones, np.ndarray):
            zones = np.array(zones)
        nzones = len(zones)

        mat = np.zeros((nzones, nzones), dtype=np.float64)
        xm = xm[zones]

        # compute contributions
        # mat of [i,j] contains contributions of zone j to zone i
        xmt = 0.
        j = -1
        i = 0
        for xmi in xm:
            j += 1
            xmt += xmi
            if xmt >= mbox:
                break
        if xmt < mbox:
            # complete mixing
            mat[:, :] = xm[np.newaxis, :] / xmt
        else:
            mat[i, i:j+1] = xm[i:j+1] / xmt

            for k in range(1, nzones):
                j0 = j
                xmt -= xm[k - 1]
                xmt0 = xmt
                while xmt <= mbox:
                    j += 1
                    xmt += xm[j]
                    if j == nzones - 1:
                        break
                if k <= j0:
                    mat[k, :j0+1] = mat[k-1, :j0+1] * (xmt0 / xmt)
                if j > j0:
                    mat[k, j0+1:j+1] = xm[j0+1:j+1] / xmt
                i = k
                if j == nzones - 1:
                    break
            # add flat tail
            mat[i+1:, :] = mat[i, :]
        mat = matrix_power(np.transpose(mat), iterations)

        result = copy.deepcopy(self)
        result.data[zones, ...] = np.tensordot(mat, result.data[zones, ...], axes = (0, 0))
        return result

    def mix_diff(self,
            diff,
            xm = None,
            zones = None,
            ):
        """
        Mix abundance dump and return new mixed object.

        'diff' is the mix width using diffusion-like approximation
        (Gaussian-like mixing with widths mbox).  In a diffusion
        apprximation diff would correstpond to sqrt(D * dt) where D is
        the diffusion coefficient and dt is the time [step].  Hence
        the unit of diff here is 'mass' (g) for conveniece.

        'diff' may also be an array of same dimension as xm.  The
        mixing magnitudes are defined on the zone interfaces similar
        to the mass coordinates zm.

        It appears similar results in the near field are obtained with
        a 'diff' value of about 2x the 'mix' 'mbox' value; in the far
        field, a value of sqrt(2) gives more similar results (on log
        scale).
        """
        if xm is None:
            try:
                xm = self.xm
            except AttributeError:
                raise Exception('Need to have zone mass defined.')
        assert xm.shape[0] == self.data.shape[0], 'require same number of zones in xm and data set'

        if zones is None:
            # exclude "phony" and surface zones
            if self.has_wind:
                zones = slice(self.i0, self.i1 - 1)
            else:
                zones = slice(self.i0, self.i1)

        # convert to list of zones
        if isinstance(zones, slice):
            zones = np.array(range(*zones.indices(len(xm))))
        if not isinstance(zones, np.ndarray):
            zones = np.array(zones)
        nzones = len(zones)

        xm = xm[zones]

        # compute forward matrix
        # mat of [i,j] contains contributions of zone j to zone i
        if np.isscalar(diff):
            diff = np.tile(diff, len(xm))
        elif len(diff) == len(self.xm):
            diff = diff[zones]

        assert len(diff) == len(xm)
        xmim = np.zeros_like(xm)
        xmip = np.zeros_like(xm)
        ip = np.arange(1,nzones)
        im = np.arange(0,nzones-1)
        ia = np.arange(0,nzones)

        xmi = -2 / xm
        xmip[im] = diff[im]**2 / (xm[im] + xm[im+1])
        xmim[ip] = xmip[im]

        mat = np.ndarray((3,nzones), dtype=np.float64)
        mat[0,ip] =  xmi[im] *  xmip[im]
        mat[1,ia] = -xmi[ia] * (xmip[ia] + xmim[ia]) + 1
        mat[2,im] =  xmi[ip] *  xmim[ip]

        result = copy.deepcopy(self)
        result.data[zones, ...] = solve_banded((1,1), mat, result.data[zones, ...])
        return result

    @CachedAttribute
    def abar(self):
        """
        Compute Abar as a function of coordinate
        """
        if not self.molfrac:
            A = 1 / isotope.ufunc_A(self.ions)
        else:
            A = np.ones(self.data.shape[1])
        xA = np.tensordot(self.data, A, axes = (1, 0))
        return 1 / xA

    # @CachedAttribute
    # def muI(self):
    #     """
    #     mean molecular weight per ion
    #     """
    #     if not self.molfrac:
    #         A = isotope.ufunc_A(self.ions)
    #     else:
    #         A = np.ones(self.data.shape[1])
    #     xA = np.tensordot(self.data, A, axes = (1, 0))
    #     return 1 / xA

    muI = abar

    @CachedAttribute
    def zbar(self):
        """
        mean charge number
        """
        Z = isotope.ufunc_Z(self.ions)
        if not self.molfrac:
            A = isotope.ufunc_A(self.ions)
            Z /= A
        else:
            A = np.ones(self.data.shape[1])
        xZ = np.tensordot(self.data, Z, axes = (1, 0))
        xA = np.tensordot(self.data, A, axes = (1, 0))
        return xZ / xA

    @CachedAttribute
    def Ye(self):
        """
        Compute Ye as a function of coordinate
        """
        Z = isotope.ufunc_Z(self.ions)
        if not self.molfrac:
            Z /= isotope.ufunc_A(self.ions)
        Ye = np.tensordot(self.data, Z, axes = (1, 0))
        return Ye

    @CachedAttribute
    def eta(self):
        """
        neutron excess eta = 1-2*Ye
        """
        return 1 - 2 * self.Ye

    @CachedAttribute
    def mu(self):
        """
        mean molecular weight
        """
        Z1 = 1 + isotope.ufunc_Z(self.ions)
        if not self.molfrac:
            Z1 /= isotope.ufunc_A(self.ions)
        xmu = np.tensordot(self.data, Z1, axes = (1, 0))
        return 1 / xmu

    @CachedAttribute
    def mue(self):
        """
        mean molecular weight per electron
        """
        return 1 / self.Ye

    @CachedAttribute
    def zm_sun(self):
        """
        mass shell interfaces in solar masses
        """
        return self.zm / physconst.XMSUN

########################################################################

class AbuSet(Logged):
    """
    Prototype Class for single abundace set like solar abundance

    Should there be a subclass or a field for normalized/partial abundances?
    - Normalization is difficult if addition occurs piecewise ...

    How about support for different formats, like [] ... ?

    *** TODO ***
    - derive using common functions with IonSet
    - replace internal variables with "_"-names,
      most prominently: iso --> _ions
                        abu --> _abu

    *** Add a way to hold non-normalized values?
    """
    def __init__(self,
                 iso = None,
                 abu = None,
                 comment = None,
                 mixture = None,
                 dat_file = None,
                 bg_file = None,
                 silent = False,
                 # todo
                 normalize = False, # only useful for fractions
                 normalized = False, # only useful for fractions
                 molfrac = False, #?
                 unit = 'massfrac', # set up proper unit system
                 # - g, mol, massfrac, molfrac, M_sun
                 sorted = False,
                 sentinel = None,
                 comp = None,
                 **kwargs):
        """
        Initialize abundance set.

        Currently allow call with
         * list of isotopes only (abu set to 0)
         * list of isotopes and list of abundances
         * keyword iso=abu
         * dictionary {'iso':abu, ...}
        Initialize from data file:
         * dat_file (KEPLER abundance data file, like solabu.dat)
         * bg_file (KEPLER BURN generator data file)

        TODO:
          * sorting implementation
          * implement molfrac
          * proper deal with normalization on/off
          * make A, Z, etc., properties rather than functions
          * can initialise from another composition for derived classes?
        """
        self.comment = stuple(comment)
        self.mixture = mixture
        self.sentinel = sentinel
        self.comp = comp
        # TODO - implement
        assert molfrac == False
        self.molfrac = molfrac
        self.unit = unit
        self.is_sorted = sorted
        if dat_file is not None:
            self._from_dat(dat_file, silent = silent)
            return
        if bg_file is not None:
            self._from_bg(bg_file, silent = silent)
            return
        if issubclass(dict, type(iso)):
            self.iso, self.abu = self._ion_abu_from_dict(**iso)
        elif iso is None:
            assert abu is None, "Need isotope name"
            self.iso = np.array([], dtype=np.object)
            self.abu = np.array([], dtype=np.float64)
        else:
            self.iso = np.array([isotope.ion(i) for i in np.atleast_1d(iso)], dtype=np.object)
            if abu is not None:
                self.abu = np.array(np.atleast_1d(abu), dtype=np.float64)
                assert len(self.abu) == len(self.iso), "Need equal number of elements."
            else:
                self.abu = np.zeros_like(self.iso, dtype=np.float64)
        if len(kwargs) > 0:
            self._append(*self._ion_abu_from_dict(**kwargs))
        if self.is_sorted:
            self.sort()
        # need to update for general use, e.g., mapped elemental abundances
        if molfrac:
            self.abu *= self.A()

    @staticmethod
    def _ion_abu_from_dict(**kwargs):
        #todo: add sorting?
        #todo: check uniqueness
        return (
            np.array([isotope.ion(i) for i in kwargs.keys()],
                     dtype=np.object),
            np.array(list(kwargs.values()),
                     dtype=np.float64))

    def _append(self, iso, abu):
        #todo: add sorting? - YES
        #todo: check uniqueness
        self.iso = np.append(self.iso, iso)
        self.abu = np.append(self.abu, abu)

    def __str__(self):
        return ("abu([" + self.unit + "] " +
                ", ".join(['{:s}: {:8G}'\
                           .format(iso.Name(),abu)
                           for iso, abu in zip(self.iso, self.abu)]) +
                ")")
    __repr__ = __str__

    def _delete(self, iso):
        """
        remove isotope
        """
        ii = np.where(self.iso == iso)[0]
        assert len(ii) == 1
        self.iso = np.delete(self.iso, ii)
        self.abu = np.delete(self.abu, ii)

    def __delitem__(self, iso):
        """
        remove isotope
        """
        self._delete(iso)

    def normalize(self, total = None):
        """
        Normalize abundances to one.

        If sum == 0 just return.
        """
        abusum = self.abu.sum()
        if abusum == 0.:
            return
        self.abu /= abusum
        if total is not None:
            self.abu *= total

    def normalized(self, total = None):
        """
        Return normalized copy of abu.
        """
        x = copy.copy(self)
        x.normalize(total = total)
        return x

    def _from_bg(self,
                 file_name,
                 silent = False):
        """
        Generate abundace set from BURN gen file.

        TODO - return tuple with mutiple mixtures if file has several

               Actually, this should be a module function, not part of
               the constructor.

        TODO - add option to show comment
        """
        self.setup_logger(silent = silent)
        self.comment = stuple(self.comment)
        xre = re.compile('[-+a-zA-Z0-9.]+')
        self.iso = np.array([],dtype=np.object)
        self.abu = np.array([],dtype=np.float64)
        with open(file_name,'r') as f:
            self.logger_file_info(f)
            for line in f:
                if line.startswith('c '):
                    self.comment += (line[2:].rstrip(),)
                elif line.startswith('m '):
                    xdata = xre.findall(line)
                    xnum = len(xdata)
                    if xnum < 2:
                        continue
                    mixture = xdata[1]
                    if self.mixture is None:
                        self.mixture = mixture
                        self.logger.info('Loading mixture "{:s}".'.format(self.mixture))
                    if self.mixture == mixture:
                        if xnum < 3:
                            continue
                        assert xnum % 2 == 0
                        xion = xdata[3::2]
                        xabu = xdata[2::2]
                        for i,a in zip(xion,xabu):
                            self._append(isotope.ion(i),np.double(a.replace('D','E').replace('d','e')))
        self.close_logger(timing='Data loaded in')


    def _from_dat(self,
                  filename,
                  silent = False):
        """
        Load abundace set from "dat" file.

        TODO - add option to show comment
        """
        self.setup_logger(silent = silent)
        self.comment = stuple(self.comment)
        xre = re.compile('[-+a-zA-Z0-9.]+')
        self.iso = np.array([],dtype=np.object)
        self.abu = np.array([],dtype=np.float64)
        with open(filename,'r') as f:
            self.logger_file_info(f)
            self.comment += ('',
                             'Generated from file "{:s}".'.format(filename),
                             'Original file comments follow:',
                             '')
            for line in f:
                if not line.startswith((';','#')):
                    xdata = xre.findall(line)
                    xnum = len(xdata)
                    if xnum == 0:
                        continue
                    if xnum == 2:
                        xion,xabu = tuple(xdata)
                    else:
                        print(line)
                        raise IOError('bad format')
                    self._append(isotope.ion(xion), np.double(xabu))
                else:
                    self.comment += (line[2:].rstrip(),)
        self.filename = filename
        message = "{:3d} isotopes loaded in".format(len(self.iso))
        self.close_logger(timing = message)

    def _new_order(self):
        # reset things that depend on order
        # this should be called every time things are added or removed.
        # try:
        #     del self.approx_map
        # except:
        #     pass
        pass

    def ionset(self):
        return IonSet(self.ions)

    def sort(self):
        """
        Sort ions.
        """
        if self.iso.size == 0:
            return
        sort = self.iso.argsort()
        self.iso = self.iso[sort]
        self.abu = self.abu[sort]
        self._new_order()

    def write_bg(self,
                 outfile,
                 net = 1,
                 mixture = None,
                 zmax = 83,
                 overwrite = False,
                 silent = False,
                 write_net = True):
        """
        Write out BURN generator.

        If outfile is file use this.
        If outfile is a filename open outfile for writing.

        We need to assert gaps around be8, b9
        (though it is not clear why KEPLER needs these gaps)

        We do need to assert, however, because ADAPNET preserves gaps
        for this reason, not to introduce gaps in new network.  Which
        also makes the work of ADAPNET easier.

        Eventually KEPLER should check for these gaps and issue errors
        if present.  Maybe it already does?
        """

        self.setup_logger(silent = silent)

        # minnet from adapnet general purpose network
        # we could load this, but for now, we just hard-code.
        # note the gaps for Be8 and B9 that appear to be requires
        # by KEPLER - though it is not clear why
        # maybe because these istopes are not in bdat
        # TODO - ask RDH to add to future versions of bdat.
        minnet = [[[ 1, 1]],
                  [[ 1, 3]],
                  [[ 3, 4]],
                  [[ 6, 7]],
                  [[ 7, 7], [ 9, 9]],
                  [[ 8, 8], [10, 11]],
                  [[11,13]],
                  [[13,15]],
                  [[14,18]],
                  [[17,19]],
                  [[19,22]],
                  [[21,23]],
                  [[23,26]],
                  [[25,27]],
                  [[27,30]],
                  [[30,31]],
                  [[31,36]]]

        version = 10000

        default_mixture = 'x'

        card_cmt = 'c'
        card_mix = 'm'
        card_grd = 'gg'
        card_ntw = 'netw'

        zmax = min(zmax+1, len(isotope.elements))
        netw = np.ndarray((zmax,2,2), dtype=np.int64)
        netw_count = np.zeros(zmax, dtype=np.int64)
        for iz,z in enumerate(minnet):
            for ia,a in enumerate(z):
                netw[iz,ia,:] = a
            netw_count[iz] = len(z)

        niso = 0
        iso = np.zeros_like(self.iso)
        abu = np.zeros_like(self.abu)
        for i,a in zip(self.iso, self.abu):
            Z = i.Z
            if Z < zmax:
                A = i.A
                n = netw_count[Z]
                if n == 0:
                    netw_count[Z] = 1
                    netw[Z,0,:] = A
                else:
                    netw[Z,0,0] = min(A, netw[Z,0,0])
                    netw[Z,n-1,1] = max(A, netw[Z,n-1,1])
                add = True
                if n  == 2:
                    if (A > netw[Z,0,1]) and (A < netw[Z,1,0]):
                        self.logger.error('{:s} inside required gap.'.format(i.name()))
                        add = False
                if add:
                    iso[niso] = i
                    abu[niso] = a
                    niso += 1
        iso = iso[:niso]
        abu = abu[:niso]

        if mixture is None:
            try:
                mixture = self.mixture
            except:
                pass
        if mixture is None:
            mixture = default_mixture

        if not isinstance(outfile, io.IOBase):
            filename = os.path.expanduser(os.path.expandvars(outfile))
            assert overwrite or not os.path.exists(filename),'file exisits'
            f = open(filename,'w')
        else:
            f = outfile

        if write_net:
            f.write(card_cmt + ' COMPUTER-GENERATED BURN GENERATOR FILE\n')
            f.write(card_cmt + ' VERSION {:s}'.format(version2human(version))+'\n')
            f.write(card_cmt + ' ' + time.asctime(time.gmtime())+' UTC\n')
            f.write(card_cmt + '\n')

            for c in self.comment:
                f.write("{:s} {:s}\n".format(card_cmt,c))
                f.write(card_cmt + '\n')
                f.write(card_cmt + ' define network (Z_max = {:d})\n'.format(zmax-1))
            for i in range(zmax):
                nc = netw_count[i]
                if nc > 0:
                    c = " ".join(["{:3d} {:3d}".format(*netw[i,j,:]) for j in range(nc)])
                    f.write("{:s} {:d} {:2s} {:s}\n".format(
                        card_ntw,
                        net,
                        isotope.elements[i],
                        c,
                        ))
        f.write(card_cmt + '\n')
        f.write(card_cmt + ' define composition '+mixture+'\n')
        for i,a in zip(iso, abu):
            f.write("{:s} {:s} {:14.8E} {:s}\n".format(card_mix,mixture,a,i.name()))
        if write_net:
            f.write(card_cmt + '\n')
            f.write(card_cmt + ' specify grid composition (homogeneous star)\n')
            f.write(card_cmt + ' NOTE: call only after all "g" cards in main generator\n')
            f.write("{:s} {:d} {:s}\n".format(card_grd,net,mixture))

        if not isinstance(outfile, io.IOBase):
            f.close()

        if write_net:
            self.close_logger(timing='BURN generator written to "{:s}" in'.format(f.name))
        else:
            self.close_logger(timing='BURN mixture {:s} written in'.format(mixture))

    def write_dat(self,
                  outfile,
                  overwrite = False,
                  silent = False):
        """
        Write out dat file.

        If file is file use this.
        If file is file name open file for write.
        """
        self.setup_logger(silent = silent)

        version = 10000

        card_cmt = ';'

        if not isinstance(outfile, io.IOBase):
            filename = os.path.expanduser(os.path.expandvars(outfile))
            assert overwrite or not os.path.exists(filename),' file exisits: '+filename
            f = open(filename,'w')
        else:
            f = outfile

        f.write(card_cmt + ' COMPUTER-GENERATED ABUNDANCE DATA FILE\n')
        f.write(card_cmt + ' VERSION {:s}'.format(version2human(version))+'\n')
        f.write(card_cmt + ' ' + time.asctime(time.gmtime())+' UTC\n')
        f.write(card_cmt + '\n')

        for c in self.comment:
            f.write("{:s} {:s}\n".format(card_cmt,c))
        f.write(card_cmt + '\n')
        f.write(card_cmt + '----------------------------------------------------------------------\n')
        for i,a in zip(self.iso, self._X()):
            f.write("{:6s} {:13.7E}\n".format(i.name(),a))
        if not isinstance(outfile, io.IOBase):
            f.close()
        self.close_logger(timing = '"{:s}" written in'.format(f.name))

    def write_compsurb(self,
                  outfile,
                  overwrite = False,
                  silent = False):
        """
        Write out data file for setting compsurf BURN values.

        If file is file use this.
        If file is file name open file for write.
        """
        self.setup_logger(silent = silent)

        version = 10000

        card_cmt = 'c'

        if not isinstance(outfile, io.IOBase):
            filename = os.path.expanduser(os.path.expandvars(outfile))
            assert overwrite or not os.path.exists(filename),' file exisits: '+filename
            f = open(filename,'w')
        else:
            f = outfile

        f.write(card_cmt + ' COMPUTER-GENERATED BURN COMPSURF LINK FILE\n')
        f.write(card_cmt + ' VERSION {:s}'.format(version2human(version))+'\n')
        f.write(card_cmt + ' ' + time.asctime(time.gmtime())+' UTC\n')
        f.write(card_cmt + '\n')

        for c in self.comment:
            f.write("{:s} {:s}\n".format(card_cmt,c))
        f.write(card_cmt + '\n')
        f.write(card_cmt + '----------------------------------------------------------------------\n')
        f.write("compsurb clear\n")
        for i,a in zip(self.iso, self._X()):
            f.write("compsurb {:13.7E} {:s}\n".format(a, i.name()))
        f.write("mapsurfb\n")
        if not isinstance(outfile, io.IOBase):
            f.close()
        self.close_logger(timing = '"{:s}" written in'.format(f.name))

    # some property routines
    # ** TODO replace with isotope.ufunc_* routines

    def A(self):
        return isotope.ufunc_A(self.iso)
    def Z(self):
        return isotope.ufunc_Z(self.iso)
    def N(self):
        return isotope.ufunc_N(self.iso)
    def E(self):
        return isotope.ufunc_E(self.iso)

    def isomer(self):
        return isotope.ufunc_isomer(self.iso)
    def isotope(self):
        return isotope.ufunc_isotope(self.iso)
    def element(self):
        return isotope.ufunc_element(self.iso)
    def isobar(self):
        return isotope.ufunc_isobar(self.iso)
    def isotone(self):
        return isotope.ufunc_isotone(self.iso)

    def element_name(self):
        return np.array([x.element_name() for x in self.iso])
    def element_symbol(self, upcase = True):
        return np.array([x.element_symbol(upcase = upcase) for x in self.iso])
    def idx(self):
        return isotope.ufunc_idx(self.iso)

    def isotone_idx(self):
        return isotope.ufunc_isotone_idx(self.iso)
    def isobar_idx(self):
        return isotope.ufunc_isobar_idx(self.iso)
    def element_idx(self):
        return isotope.ufunc_element_idx(self.iso)
    def isotope_idx(self):
        return isotope.ufunc_isotope_idx(self.iso)
    def isomer_idx(self):
        return isotope.ufunc_isomer_idx(self.iso)

    def isotones_idx(self):
        return np.unique(self.isotone_idx())
    def isobars_idx(self):
        return np.unique(self.isobar_idx())
    def elements_idx(self):
        return np.unique(self.element_idx())
    def isotopes_idx(self):
        return np.unique(self.isotope_idx())
    def isomers_idx(self):
        return np.unique(self.isomer_idx())

    def get_isotopes(self, ion):
        ion = isotope.ion(ion)
        ii = np.where(isotope.ufunc_Z(self.iso) == ion.Z)[0]
        return self.iso[ii]
    def get_isotones(self, ion):
        ion = isotope.ion(ion)
        ii = np.where(isotope.ufunc_N(self.iso) == ion.N)[0]
        return self.iso[ii]
    def get_isobars(self, ion):
        ion = isotope.ion(ion)
        ii = np.where(isotope.ufunc_A(self.iso) == ion.A)[0]
        return self.iso[ii]
    def get_isomers(self, ion):
        ion = isotope.ion(ion)
        ii = np.where(np.logical_and(
            isotope.ufunc_Z(self.iso) == ion.Z,
            isotope.ufunc_N(self.iso) == ion.N)
            )[0]
        return self.iso[ii]

    def XYZ(self):
        """
        Return 'astronomical' X, Y, Z of composition by mass fraction
        """
        x = sum(self.X()[self.Z() <= 1])
        y = sum(self.X()[self.Z() == 2])
        z = sum(self.X()[self.Z() >= 3])
        return np.array([x, y, z])

    def xyz(self):
        """
        Return 'astronomical' X, Y, Z of composition by mol fraction
        """
        x = sum(self.Y()[self.Z() <= 1])
        y = sum(self.Y()[self.Z() == 2])
        z = sum(self.Y()[self.Z() >= 3])
        return np.array([x, y, z])

    def metallicity(self):
        """
        Return 'metallicity' Z of composition.
        """
        z = sum(self.X()[self.Z() >= 3])
        return z

    def Ye(self):
        """
        Return electron to baryon ratio.
        """
        Ye = np.sum(self.Z() * self.Y())
        return Ye

    def mue(self):
        """
        Return mean molecular weight per electron of composition.
        """
        return 1 / self.Ye()

    def eta(self):
        """
        Return neutron excess of mixture.
        """
        return 1 - 2 * self.Ye()

    def mu(self):
        """
        Return mean molecular weight of composition.
        """
        xmu = np.sum((self.Z() + 1) * self.Y())
        return 1 / xmu

    # here to add some general routines that return data for any Ion type
    # (isotope, isotone, isobar, element)
    def _get(self, selection = None, data = 'X', invalid = np.nan, exception = False):
        """
        Return data for selected (isotope, isotone, isobar, element).

        If nothing is specified, all isotopes will be returned.
        Default is mass fraction 'X'.
        """

        def output(x):
            nonlocal shape, exception, selection
            if exception:
                index, = np.where(np.isnan(x))
                print(index)
                if len(index) > 0:
                    raise AttributeError('Isotope(s) not found: ' + ', '.join([str(ion) for ion in selection[index]]))
            if len(shape) == 0:
                x = x[0]
            elif len(shape) != 1:
                x = x.reshape(shape)
            return x

        # the following should go into class definition or init routine

        # setup
        uis_func = [
            isotope.ufunc_is_isomer,
            isotope.ufunc_is_isotope,
            isotope.ufunc_is_element,
            isotope.ufunc_is_isobar,
            isotope.ufunc_is_isotone,
            ]

        # we would only have to switch out this one definition to get them all...
        if data == 'X':
            val_func = [
                self.isomers_X,
                self.isotopes_X,
                self.elements_X,
                self.isobars_X,
                self.isotones_X,
                ]
        elif data == 'Y':
            val_func = [
                self.isomers_Y,
                self.isotopes_Y,
                self.elements_Y,
                self.isobars_Y,
                self.isotones_Y,
                ]
        elif data == 'log_eps':
            val_func = [
                self.isomers_log_eps,
                self.isotopes_log_eps,
                self.elements_log_eps,
                self.isobars_log_eps,
                self.isotones_log_eps,
                ]
        elif data == 'Abu':
            val_func = [
                self.isomers_Abu,
                self.isotopes_Abu,
                self.elements_Abu,
                self.isobars_Abu,
                self.isotones_Abu,
                ]
        else:
            raise ValueError("Invalid Data request: '{}'.".format(data))

        # selection setup
        if selection is None:
            selection = self.iso
        shape = np.shape(selection)
        if shape == ():
            selection = np.array([selection])
        if not isinstance(selection, np.ndarray):
            selection = np.array(selection)
        if len(shape) > 1:
            selection = selection.reshape(-1)
        # use factory function instead of Ion
        selection = np.array([isotope.ion(ix) for ix in selection])

        selections = np.ndarray([len(uis_func), len(selection)], dtype = np.bool)

        # pure cases
        # here we test the most likely cases first ...
        for i,(uis,val) in enumerate(zip(uis_func, val_func)):
            selections[i, :] = uis(selection)
            if np.all(selections[i, :]):
                return output(val(selection))

        # now mixed cases
        if exception:
            invalid = np.nan
        x = np.tile(np.array(invalid, dtype=np.float64), len(selection))
        for i,val in enumerate(val_func):
            if np.any(selections[i,:]):
                index, = np.nonzero(selections[i,:])
                x[index] =  val(selection[index])

        return output(x)

    def _X(self):
        return self.abu

    def _Y(self):
        return self.abu / self.A()

    def ppm(self, selection = None):
        """
        Return ppm (10,000 ppm = 1 mass%) for selected(isomer, isotope, isotone, isobar, element).

        If nothing is specified, all isotopes will be returned.
        """
        return self.X(selection) * 1e6

    def X(self, selection = None):
        """
        Return mass fraction for selected (isomer, isotope, isotone, isobar, element).

        If nothing is specified, all isotopes will be returned.
        """
        return self._get(selection, data = 'X')

    def Y(self, selection = None):
        """
        Return mole/g for selected (isomer, isotope, isotone, isobar, element).

        If nothing is specified, all isotopes will be returned.
        """
        return self._get(selection, data = 'Y')

    def log_eps(self, selection = None):
        """
        Return log(eps) for selected (isomer, isotope, isotone, isobar, element).

        If nothing is specified, all isotopes will be returned.
        """
        return self._get(selection, data = 'log_eps')

    def Abu(self, selection = None):
        """
        Return A with A(Si) = 10**6 for selected (isomer, isotope, isotone, isobar, element).

        If nothing is specified, all isotopes will be returned.
        """
        return self._get(selection, data = 'Abu')

    # add functions to compute [] , \delta - those need to be in SolAbu


    # generic selection routine
    def _selection(self,
                   x,
                   idx,
                   selection = None,
                   missing = np.nan,
                   check = None,
                   convert = None,
                   ):
        if selection is None:
            return x
        # we should really use IonList objects
        if np.shape(selection) == ():
            selection = np.array([selection])
        selection = isotope.ufunc_ion(selection)
        if not np.all(check(selection)):
            raise KeyError("Selection contained invalid entries.")
        selection = convert(selection)
        y = np.ndarray(selection.shape, dtype = x.dtype)
        ii = np.argsort(idx)
        sel = np.minimum(np.searchsorted(idx[ii], selection), len(ii) - 1)
        jj = idx[ii][sel] == selection
        y[jj] = x[ii][sel][jj]
        y[~jj] = missing
        return y


    def _generic_selection(self,
                           data = None,
                           selection = None,
                           mode = None,
                           missing = np.nan,
                           return_selection = False,
                           **kwargs):
        if data is None:
            data = self.abu
        elif data == 'X':
            data = self._X()
        elif data == 'Y':
            data = self._Y()
        if mode is None:
            convert = isotope.ufunc_idx
            check = isotope.ufunc_is_ion
            keys = self.idx()
        else:
            keys = self.__getattribute__(mode + '_idx')(),
            check = isotope.__dict__['ufunc_is_' + mode]
            convert = isotope.__dict__['ufunc_' + mode + '_idx']
        x, idx = project(
            data,
            keys,
            return_values = True,
            )
        if selection is None:
            result = x
        else:
            result = self._selection(
                x,
                idx,
                selection = selection,
                check = check,
                convert = convert,
                missing = missing,
                **kwargs)
        if return_selection:
            if selection is None:
                selection = isotope.ufunc_ion_from_idx(idx)
            return result, selection
        return result

    # ------------------------------
    # isomer routines
    # ------------------------------
    def _isomers_selection(self, **kwargs):
        return self._generic_selection(
            mode = 'isomer',
            **kwargs)

    def isomers_X(self, selection = None):
        """
        Return isomer mass fractions.
        """
        return self._isomers_selection(
            data = 'X',
            selection = selection)

    def isomers_Y(self, selection = None):
        """
        Return isomer number fractions (mol per gram).
        """
        return self._isomers_selection(
            data = 'Y',
            selection = selection)

    def isomers_log_eps(self, selection = None):
        """
        Return isomer log(eps).

        log10(# of atoms relative to H) + 12
        """
        h = self.Z() == 1
        A = self.A()
        x = sum(self._X()[h] / A[h])
        y = self.isomers_Y(selection)
        return np.log10(y / x) + 12

    def isomers_set(self, selection = None):
        """
        Return set of projected isomers.
        """
        abu, iso = self._isomers_selection(
            data = 'X',
            selection = selection,
            return_selection = True)
        return AbuSet(abu = abu,
                      iso = iso,
                      comment = self.comment)
        # we may want to add something signifying isomers ...

    def isomers_log_eps(self, selection = None):
        """
        Return isomere log(eps).

        log10(# of atoms relative to H) + 12
        """
        h = self.Z() == 1
        A = self.A()
        x = sum(self._X()[h] / A[h])
        y = self.isomers_Y(selection)
        return np.log10(y / x) + 12

    def isomers_Abu(self, selection = None):
        """
        Return isomere Abundance A(Si) = 10**6.

        number of atoms relative to Si
        """
        si = self.Z() == 14
        A = self.A()
        x = sum(self._X()[si] / A[si])
        y = self.isomers_Y(selection)
        return y / x * 1e6

    def isomers(self):
        """
        Return list of all isomers in abundance set.
        """
        idx = np.unique(self.isomer_idx())
        ii = np.where(idx != isotope.Ion.VOID_IDX)
        return isotope.ufunc_ion_from_idx(idx[ii])

    # ------------------------------
    # isotope routines
    # ------------------------------
    def _isotopes_selection(self, **kwargs):
        return self._generic_selection(mode = 'isotope', **kwargs)

    def isotopes_X(self, selection = None):
        """
        Return isotope mass fractions.
        """
        return self._isotopes_selection(
            data = 'X',
            selection = selection)

    def isotopes_Y(self, selection = None):
        """
        Return isotope number fractions (mol per gram).
        """
        return self._isotopes_selection(
            data = 'Y',
            selection = selection)

    def isotopes_log_eps(self, selection = None):
        """
        Return isotope log(eps).

        log10(# of atoms relative to H) + 12
        """
        h = self.Z() == 1
        A = self.A()
        x = sum(self._X()[h] / A[h])
        y = self.isotopes_Y(selection)
        return np.log10(y / x) + 12


    def isotopes_Abu(self, selection = None):
        """
        Return isotope Abundance A(Si) = 10**6.

        number of atoms relative to Si
        """
        si = self.Z() == 14
        A = self.A()
        x = sum(self._X()[si] / A[si])
        y = self.isotopes_Y(selection)
        return y / x * 1e6


    def isotopes(self):
        """
        Return list of all isotopes in abundance set.
        """
        idx = np.unique(self.isotope_idx())
        ii = np.where(idx != isotope.Ion.VOID_IDX)
        return isotope.ufunc_ion_from_idx(idx[ii])


    def isotopes_set(self, selection = None):
        """
        Return set of projected isotopes.
        """
        abu, iso = self._isotopes_selection(
            data = 'X',
            selection = selection,
            return_selection = True)
        return AbuSet(abu = abu,
                      iso = iso,
                      comment = self.comment)
        # we may want to add something signifying isotopes ...


    # ------------------------------
    # element routines
    # ------------------------------
    def _elements_selection(self, **kwargs):
        return self._generic_selection(mode = 'element', **kwargs)

    def elements(self):
        """
        Return list of all elements in abundance set.
        """
        idx = np.unique(self.element_idx())
        ii = np.where(idx != isotope.Ion.VOID_IDX)
        return isotope.ufunc_ion_from_idx(idx[ii])

    def elements_Z(self):
        """
        Return charge number of all elements.
        """
        return np.unique(self.Z())

    def elements_X(self, selection = None):
        """
        Return elemental mass fractions.
        """
        return self._elements_selection(
            data = 'X',
            selection = selection)

    def elements_Y(self, selection = None):
        """
        Return elemental number fractions (mol per gram).
        """
        return self._elements_selection(
            data = 'Y',
            selection = selection)

    def elements_log_eps(self, selection = None):
        """
        Return elemental log(eps).

        log10(# of atoms relative to H) + 12
        """
        h = self.Z() == 1
        A = self.A()
        x = sum(self._X()[h] / A[h])
        y = self.elements_Y(selection)
        return np.log10(y / x) + 12

    def elements_Abu(self, selection = None):
        """
        Return element abundance A(Si) = 10**6.

        number of atoms relative to Si
        """
        si = self.Z() == 14
        A = self.A()
        x = sum(self._X()[si] / A[si])
        y = self.elements_Y(selection)
        return y / x * 1e6

    def elements_name(self):
        """
        Return element name of each isotope.
        """
        return np.array([Elements[x] for x in self.elements_Z()])

    def elements_set(self, selection = None):
        """
        Return set of projected elements.
        """
        abu, iso = self._elements_selection(
            data = 'X',
            selection = selection,
            return_selection = True)
        return AbuSet(abu = abu,
                      iso = iso,
                      comment = self.comment)
        # we may want to add something signifying elements ...

    # ------------------------------
    # isobar routines
    # ------------------------------
    def _isobars_selection(self, **kwargs):
        return self._generic_selection(mode = 'isobar', **kwargs)

    def isobars(self):
        """
        Return list of all elements in abundance set.
        """
        idx = np.unique(self.isobar_idx())
        ii = np.where(idx != isotope.Ion.VOID_IDX)
        return isotope.ufunc_ion_from_idx(idx[ii])

    def isobars_A(self):
        """
        Return mass number of all isobars.
        """
        return np.unique(self.A())

    def isobars_X(self, selection=None):
        """
        Return isobar mass fractions.
        """
        return self._isobars_selection(
            data = 'X',
            selection = selection)

    def isobars_Y(self, selection=None):
        """
        Return isobar number fractions (mol per gram).
        """
        return self._isobars_selection(
            data = 'Y',
            selection = selection)

    def isobars_log_eps(self, selection = None):
        """
        Return isobar log(eps).

        log10(# of atoms relative to H) + 12
        """
        h = self.Z() == 1
        A = self.A()
        x = sum(self._X()[h] / A[h])
        y = self.isobars_Y(selection)
        return np.log10(y / x) + 12

    def isobars_Abu(self, selection = None):
        """
        Return isobar abundance A(Si) = 10**6.

        number of atoms relative to Si
        """
        si = self.Z() == 14
        A = self.A()
        x = sum(self._X()[si] / A[si])
        y = self.isobars_Y(selection)
        return y / x * 1e6

    def isobars_set(self, selection = None):
        """
        Return set of projected isobars.
        """
        abu, iso = self._isobars_selection(
            data = 'X',
            selection = selection,
            return_selection = True)
        return AbuSet(abu = abu,
                      iso = iso,
                      comment = self.comment)
        # we may want to add something signifying isobars ...

    # ------------------------------
    # isotone routines
    # ------------------------------
    def _isotones_selection(self, **kwargs):
        return self._generic_selection(mode = 'isotone', **kwargs)

    def isotones(self):
        """
        Return list of all isotones in abundance set.
        """
        idx = np.unique(self.isotone_idx())
        ii = np.where(idx != isotope.Ion.VOID_IDX)
        return isotope.ufunc_ion_from_idx(idx[ii])

    def isotones_N(self):
        """
        Return neutron number of all isotones.
        """
        return np.unique(self.N())

    def isotones_X(self, selection=None):
        """
        Return isotone mass fractions.
        """
        return self._isotones_selection(
            data = 'X',
            selection = selection)

    def isotones_Y(self, selection=None):
        """
        Return isotone number fractions (mol per gram).
        """
        return self._isotones_selection(
            data = 'Y',
            selection = selection)

    def isotones_log_eps(self, selection = None):
        """
        Return isotone log(eps).

        log10(# of atoms relative to H) + 12
        """
        h = self.Z() == 1
        A = self.A()
        x = sum(self._X()[h] / A[h])
        y = self.isotones_Y(selection)
        return np.log10(y / x) + 12

    def isotones_Abu(self, selection = None):
        """
        Return isotone abundance A(Si) = 10**6.

        number of atoms relative to Si
        """
        si = self.Z() == 14
        A = self.A()
        x = sum(self._X()[si] / A[si])
        y = self.isotones_Y(selection)
        return y / x * 1e6

    def isotones_set(self, selection = None):
        """
        Return set of projected isotones.
        """
        abu, iso = self._isotones_selection(
            data = 'X',
            selection = selection,
            return_selection = True)
        return AbuSet(abu = abu,
                      iso = iso,
                      comment = self.comment)
        # we may want to add something signifying isotones ...

    # general access interfaces
    def __getitem__(self, index):
        try:
            return self._get(index)
        except:
            raise AttributeError('Isotope not found.')

    def __setitem__(self, index, item):
        # TODO add isotope if  not in list?
        # maybe add parameter allow_new
        try:
            # this does not work for numpy
            if isinstance(index, str):
                index = isotope.ion(index)
            if isinstance(index, isotope.Ion):
                index = np.where(self.iso == index)
            if not is_numlike(item):
                raise AttributeError('Abundance needs to be numlike.')
            self.abu[index] = item
        except:
            raise AttributeError('Isotope not found.')

    def __len__(self):
        """
        Return number of isotopes.
        """
        return len(self.iso)

    def __getattr__(self, attr):
        try:
            x = self._get(attr, invalid = np.nan)
        except AttributeError:
            pass
        else:
            if not np.isnan(x):
                return x
        return super().__getattribute__(attr)

    def __setattr__(self, attr, value):
        if attr not in self.__dict__:
            if 'iso' in self.__dict__:
                x = np.where(self.iso == attr)[0]
                if len(x) == 1:
                    i = x[0]
                    fac = 1 - value
                    self.abu[i] = 0.
                    self._normalize(fac) # not implemented
                    self.abu[i] = value
                    return
        super().__setattr__(attr,value)

    def __iter__(self):
        i=0
        while i < self.__len__():
            yield (self.iso[i], self.abu[i])
            i+=1

    def __call__(self, *args, **kwargs):
        """
        TODO:
        1) add new isotopes to abundance pattern
        2) renormailzie
        3) check value is in Range 0<v<1
        """
        for k in args:
            x = np.where(self.iso == k)[0]
            if len(x) == 1:
                return self.abu[x[0]]
            else:
                raise AttributeError('Isotope not in list')
        for k,v in kwargs.items():
            if not 0 <= v <= 1:
                raise ValueError('Abundance for '+k+' is out of range.')
            x=np.where(self.iso == k)[0]
            if len(x) == 1:
                self.abu[x[0]] = v
            else:
                # TODO - add isotope instead (if valid)
                raise AttributeError('Isotope not in list')


    def __contains__(self, other):
        """
        Determine whether AbuSet contains other (iso) or AbuSet (vector).
        """
        if isinstance(other, AbuSet):
            return self.contains(other).all()
        x, = np.where(self.iso == other)
        if len(x) == 1:
            return True
        return False

    def contains(self, iso):
        """
        Determine whether AbuSet contains iso (scalar) or which isotopes in iso are present (vector).

        If scalar argument: return True/False.
        If vector argument return array with result for each element.
        """
        if np.isscalar(iso):
            return iso in self
        else:
            if isinstance(iso, AbuSet):
                iso = iso.iso
            return self._in1d(iso, self.iso)

    @staticmethod
    def _in1d(ar1, ar2):

        """
        replace non-working in1d that uses mergesort

        we convert to index and then call the original in1d
        ... which is sort of what would have happened anyway ...
        """
        ar1 = np.array([i.index() for i in ar1])
        ar2 = np.array([i.index() for i in ar2])
        return np.in1d(ar1, ar2)
        # ar = np.concatenate( (ar1, ar2) )
        # order = ar.argsort()
        # sar = ar[order]
        # equal_adj = (sar[1:] == sar[:-1])
        # flag = np.concatenate( (equal_adj, [False] ) )
        # indx = order.argsort()[:len( ar1 )]
        # return flag[indx]

    def ions(self):
        """
        return IonList
        """
        return IonList(self.iso)

    # ------------------------------
    # arithmetic operations
    # ------------------------------
    _zero_missing = None

    @staticmethod
    def _return_matching(iso1, abu1, iso2, abu2,
                         missing = None,
                         missing1 = None,
                         missing2 = None):
        """
        if missing is None: return intersection of isotopes, otherwise union

        if only one missing value is provided, use other set as basis
        """
        if missing1 is None:
            missing1 = missing
        if missing2 is None:
            missing2 = missing
        idx1 = isotope.ufunc_idx(iso1)
        idx2 = isotope.ufunc_idx(iso2)
        if missing1 is None and missing2 is None:
            idx = np.sort(np.intersect1d(idx1, idx2))
            ii1 = np.argsort(idx1)
            ii2 = np.argsort(idx2)
            abu1 = abu1[ii1]
            abu2 = abu2[ii2]
            idx1 = idx1[ii1]
            idx2 = idx2[ii2]
            ii1 = np.searchsorted(idx1, idx)
            ii2 = np.searchsorted(idx2, idx)
            abu_1 = abu1[ii1]
            abu_2 = abu2[ii2]
        elif missing1 is None and missing2 is not None:
            ii1 = np.argsort(idx1)
            abu_1 = abu1[ii1]
            idx = idx1[ii1]
            abu_2 = np.ndarray(len(idx), dtype = np.float64)
            abu_2.fill(missing2)
            sel = np.minimum(np.searchsorted(idx, idx2), len(idx) - 1)
            jj = idx[sel] == idx2
            abu_2[sel[jj]] = abu2[jj]
        elif missing1 is not None and missing2 is None:
            ii2 = np.argsort(idx2)
            abu_2 = abu2[ii2]
            idx = idx2[ii2]
            abu_1 = np.ndarray(len(idx), dtype = np.float64)
            abu_1.fill(missing1)
            sel = np.minimum(np.searchsorted(idx, idx1), len(idx) - 1)
            jj = idx[sel] == idx1
            abu_1[sel[jj]] = abu1[jj]
        else:
            idx = np.sort(np.union1d(idx1, idx2))
            abu_1 = np.ndarray(len(idx), dtype = np.float64)
            abu_2 = np.ndarray(len(idx), dtype = np.float64)
            abu_1.fill(missing1)
            abu_2.fill(missing2)
            ii1 = np.searchsorted(idx, idx1)
            ii2 = np.searchsorted(idx, idx2)
            abu_1[ii1] = abu1
            abu_2[ii2] = abu2
        iso = isotope.ufunc_ion_from_idx(idx)
        return iso, abu_1, abu_2


    def __truediv__(self, other):
        """
        return abundance ratio

        TODO - check for X or Y
        """
        if isinstance(other, AbuSet):
            iso, abu1, abu2 = self._return_matching(
                self.iso,
                self.X(),
                other.iso,
                other.X(),
                missing1 = self._zero_missing)
            new = abu1 / abu2
            comment = other.comment
        else:
            try:
                new = self.abu / other
            except:
                return NotImplemented
            iso = self.iso
            comment = str(other)
        return AbuSet(
            iso,
            new,
            comment = ' / '.join(stuple(self.comment, comment))
            )

    def __rtruediv__(self, other):
        """
        return abundance ratio

        TODO - check for X or Y
        """
        if isinstance(other, AbuSet):
            iso, abu1, abu2 = self._return_matching(
                self.iso,
                self.X(),
                other.iso,
                other.X(),
                missing2 = self._zero_missing)
            new = abu2 / abu1
            comment = other.comment
        else:
            try:
                new = other / self.abu
            except:
                return NotImplemented
            comment = str(other)
            iso = self.iso
        return AbuSet(
            iso,
            new,
            comment = ' / '.join(stuple(comment, self.comment))
            )

    def __floordiv__(self, other):
        """
        abundace ratio, try working subset match

        TODO - check for X or Y
        """
        if not isinstance(other, AbuSet):
            return NotImplemented
        iso = []
        new = []
        mis = []
        for i in other.iso:
            a = self[i]
            b = other[i]
            if not np.isnan(a) and b > 0:
                iso.append(i)
                new.append(a / b)
            else:
                mis.append(i)
                break
        if len(mis) > 0:
            iso = []
            new = []
            mis = []
            for i in self.iso:
                a = self[i]
                b = other[i]
                if not np.isnan(a) and b > 0:
                    iso.append(i)
                    new.append(a / b)
                else:
                    mis.append(i)
                    break
        if len(mis) > 0:
            return NotImplemented
        return AbuSet(
            iso,
            new,
            comment = ' // '.join(stuple(self.comment, other.comment))
            )


    def __add__(self, other):
        """
        return sum of abundances

        TODO - check for X or Y
        """
        if not isinstance(other, AbuSet):
            raise NotImplementedError()
        iso, abu1, abu2 = self._return_matching(
            self.iso,
            self.X(),
            other.iso,
            other.X(),
            missing = 0)
        new = abu1 + abu2
        comment = other.comment
        return AbuSet(
            iso,
            new,
            comment = ' + '.join(stuple(self.comment, comment))
            )

    def __sub__(self, other):
        """
        return difference of abundances

        TODO - check for X or Y
        """
        if not isinstance(other, AbuSet):
            return NotImplemented
        iso, abu1, abu2 = self._return_matching(
            self.iso,
            self.X(),
            other.iso,
            other.X(),
            missing = 0)
        new = abu1 - abu2
        comment = other.comment
        return AbuSet(
            iso,
            new,
            comment = ' - '.join(stuple(self.comment, comment))
            )

    def __mul__(self, other):
        """
        return product of abundances

        TODO - check for X or Y
        """
        if isinstance(other, AbuSet):
            iso, abu1, abu2 = self._return_matching(
                self.iso,
                self.X(),
                other.iso,
                other.X(),
                missing = self._zero_missing)
            new = abu1 * abu2
            comment = other.comment
        else:
            try:
                new = self.abu * other
            except:
                return NotImplemented
            comment = str(other)
            iso = self.iso
        return AbuSet(
            iso,
            new,
            comment = ' * '.join(stuple(self.comment, comment))
            )
    __rmul__ = __mul__


    def __pow__(self, other):
        """
        return power of abundances

        TODO - check for X or Y
        """
        if isinstance(other, AbuSet):
            iso, abu1, abu2 = self._return_matching(
                self.iso,
                self.X(),
                other.iso,
                other.X(),
                missing1 = self._zero_missing)
            new = abu2 ** abu1
            comment = other.comment
        else:
            try:
                new = self.abu ** other
            except:
                return NotImplemented
            comment = str(other)
            iso = self.iso
        return AbuSet(
            iso,
            new,
            comment = ' ** '.join(stuple(self.comment, comment))
            )

    def __neg__(self):
        """
        return negative of abundance
        """
        new = - self.abu
        return AbuSet(
            self.iso,
            new,
            comment = ' - '.join(stuple('', self.comment))
            )

    def index(self, iso):
        """
        Return isotope index in iso/abu array

        accepts scalar, np.ndarray, list, tuple, set
        """
        if isinstance(iso, str) or not isinstance(iso, (Iterable, np.ndarray)):
            return np.where(self.iso == iso)[0][0]
        elif isinstance(iso, np.ndarray):
            ii = np.ndarray(iso.shape, dtype=np.int64)
            for i,k in enumerate(iso.flat):
                ii[i] = np.where(self.iso == k)[0][0]
            return ii
        elif isinstance(iso, Iterable):
            ii = []
            for k in iso:
                if isinstance(k, str):
                    i = np.where(self.iso == k)[0][0]
                else:
                    i = self.index(k)
                ii += [i]
            if not isinstance(iso, list):
                ii = type(iso)(ii)
            return ii
        else:
            raise AttributeError('Argument type not supported.')

    def kepler_write_compsurb(self, filename):
        """
        Write link file to set KEPLER surface BURN abundances.
        """
        with open(filename, 'wt') as f:
            f.write('compsurb clear\n')
            for i in self.ions():
                x = self.X(i)
                if x > 0:
                    f.write('compsurb {abu:12.5e} {ion:<5s}\n'.format(ion=i.name(), abu=x))
            f.write('compsurb show\n')
