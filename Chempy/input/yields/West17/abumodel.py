"""
provide base class for abunadnace models like scaled solar
"""

import numpy as np
from abc import abstractmethod

from isotope import ion, Ion, ufunc_Z
from abuset import AbuSet, IonList

class AbuModel(AbuSet):
    """
    Base class for abudance models

    Provides different ways to specify metallicity.

    intreface through user-defined abstract methods to provide
    abundances by mass, mol, or number fraction

        _abu_molfrac_raw(self, x)

        _abu_massfrac_raw(self, x)

        _abu_numfrac_raw(self, x):

    (only one should be defined) where x is some metallicity-like parameter
    """

    ions = IonList([
        'h1', 'h2', 'he3', 'he4', 'li6', 'li7', 'be9', 'b10', 'b11',
        'c12', 'c13', 'n14', 'n15', 'o16', 'o17', 'o18', 'f19', 'ne20',
        'ne21', 'ne22', 'na23', 'mg24', 'mg25', 'mg26', 'al27', 'si28',
        'si29', 'si30', 'p31', 's32', 's33', 's34', 's36', 'cl35', 'cl37',
        'ar36', 'ar38', 'ar40', 'k39', 'k40', 'k41', 'ca40', 'ca42', 'ca43',
        'ca44', 'ca46', 'ca48', 'sc45', 'ti46', 'ti47', 'ti48', 'ti49',
        'ti50', 'v50', 'v51', 'cr50', 'cr52', 'cr53', 'cr54', 'mn55', 'fe54',
        'fe56', 'fe57', 'fe58', 'co59', 'ni58', 'ni60', 'ni61', 'ni62',
        'ni64', 'cu63', 'cu65', 'zn64', 'zn66', 'zn67', 'zn68', 'zn70',
        'ga69', 'ga71', 'ge70', 'ge72', 'ge73', 'ge74', 'ge76', 'as75',
        'se74', 'se76', 'se77', 'se78', 'se80', 'se82', 'br79', 'br81',
        'kr78', 'kr80', 'kr82', 'kr83', 'kr84', 'kr86', 'rb85', 'rb87',
        'sr84', 'sr86', 'sr87', 'sr88', 'y89', 'zr90', 'zr91', 'zr92', 'zr94',
        'zr96', 'nb93', 'mo92', 'mo94', 'mo95', 'mo96', 'mo97', 'mo98',
        'mo100', 'ru96', 'ru98', 'ru99', 'ru100', 'ru101', 'ru102', 'ru104',
        'rh103', 'pd102', 'pd104', 'pd105', 'pd106', 'pd108', 'pd110',
        'ag107', 'ag109', 'cd106', 'cd108', 'cd110', 'cd111', 'cd112',
        'cd113', 'cd114', 'cd116', 'in113', 'in115', 'sn112', 'sn114',
        'sn115', 'sn116', 'sn117', 'sn118', 'sn119', 'sn120', 'sn122',
        'sn124', 'sb121', 'sb123', 'te120', 'te122', 'te123', 'te124',
        'te125', 'te126', 'te128', 'te130', 'i127', 'xe124', 'xe126', 'xe128',
        'xe129', 'xe130', 'xe131', 'xe132', 'xe134', 'xe136', 'cs133',
        'ba130', 'ba132', 'ba134', 'ba135', 'ba136', 'ba137', 'ba138',
        'la138', 'la139', 'ce136', 'ce138', 'ce140', 'ce142', 'pr141',
        'nd142', 'nd143', 'nd144', 'nd145', 'nd146', 'nd148', 'nd150',
        'sm144', 'sm147', 'sm148', 'sm149', 'sm150', 'sm152', 'sm154',
        'eu151', 'eu153', 'gd152', 'gd154', 'gd155', 'gd156', 'gd157',
        'gd158', 'gd160', 'tb159', 'dy156', 'dy158', 'dy160', 'dy161',
        'dy162', 'dy163', 'dy164', 'ho165', 'er162', 'er164', 'er166',
        'er167', 'er168', 'er170', 'tm169', 'yb168', 'yb170', 'yb171',
        'yb172', 'yb173', 'yb174', 'yb176', 'lu175', 'lu176', 'hf174',
        'hf176', 'hf177', 'hf178', 'hf179', 'hf180', 'ta180', 'ta181', 'w180',
        'w182', 'w183', 'w184', 'w186', 're185', 're187', 'os184', 'os186',
        'os187', 'os188', 'os189', 'os190', 'os192', 'ir191', 'ir193',
        'pt190', 'pt192', 'pt194', 'pt195', 'pt196', 'pt198', 'au197',
        'hg196', 'hg198', 'hg199', 'hg200', 'hg201', 'hg202', 'hg204',
        'tl203', 'tl205', 'pb204', 'pb206', 'pb207', 'pb208', 'bi209',
        'th232', 'u234', 'u235', 'u238',
        ])

    masses = np.array([
        1.007825, 2.014102, 3.016029, 4.002603, 6.015122, 7.016004,
        9.012182, 10.012937, 11.009305, 12.000000, 13.003355, 14.003074,
        15.000109, 15.994915, 16.999132, 17.999160, 18.998403, 19.992440,
        20.993847, 21.991386, 22.989770, 23.985042, 24.985837, 25.982593,
        26.981538, 27.976927, 28.976495, 29.973770, 30.973762, 31.972071,
        32.971458, 33.967867, 35.967081, 34.968853, 36.965903, 35.967546,
        37.962732, 39.962383, 38.963707, 39.963999, 40.961826, 39.962591,
        41.958618, 42.958767, 43.955481, 45.953693, 47.952534, 44.955910,
        45.952629, 46.951764, 47.947947, 48.947871, 49.944792, 49.947163,
        50.943964, 49.946050, 51.940512, 52.940654, 53.938885, 54.938050,
        53.939615, 55.934942, 56.935399, 57.933280, 58.933200, 57.935348,
        59.930791, 60.931060, 61.928349, 63.927970, 62.929601, 64.927794,
        63.929147, 65.926037, 66.927131, 67.924848, 69.925325, 68.925581,
        70.924705, 69.924250, 71.922076, 72.923459, 73.921178, 75.921403,
        74.921596, 73.922477, 75.919214, 76.919915, 77.917310, 79.916522,
        81.916700, 78.918338, 80.916291, 77.920386, 79.916378, 81.913485,
        82.914136, 83.911507, 85.910610, 84.911789, 86.909183, 83.913425,
        85.909262, 86.908879, 87.905614, 88.905848, 89.904704, 90.905645,
        91.905040, 93.906316, 95.908276, 92.906378, 91.906810, 93.905088,
        94.905841, 95.904679, 96.906021, 97.905408, 99.907477, 95.907598,
        97.905287, 98.905939, 99.904220, 100.905582, 101.904350, 103.905430,
        102.905504, 101.905608, 103.904035, 104.905084, 105.903483,
        107.903894, 109.905152, 106.905093, 108.904756, 105.906458,
        107.904183, 109.903006, 110.904182, 111.902757, 112.904401,
        113.903358, 115.904755, 112.904061, 114.903878, 111.904821,
        113.902782, 114.903346, 115.901744, 116.902954, 117.901606,
        118.903309, 119.902197, 121.903440, 123.905275, 120.903818,
        122.904216, 119.904020, 121.903047, 122.904273, 123.902819,
        124.904425, 125.903306, 127.904461, 129.906223, 126.904468,
        123.905896, 125.904269, 127.903530, 128.904779, 129.903508,
        130.905082, 131.904154, 133.905395, 135.907220, 132.905447,
        129.906310, 131.905056, 133.904503, 134.905683, 135.904570,
        136.905821, 137.905241, 137.907107, 138.906348, 135.907144,
        137.905986, 139.905434, 141.909240, 140.907648, 141.907719,
        142.909810, 143.910083, 145.913112, 147.916889, 149.920887,
        144.912744, 143.911995, 146.914893, 147.914818, 148.917180,
        149.917271, 151.919728, 153.922205, 150.919846, 152.921226,
        151.919788, 153.920862, 154.922619, 155.922120, 156.923957,
        157.924101, 159.927051, 158.925343, 155.924278, 157.924405,
        159.925194, 160.926930, 161.926795, 162.928728, 163.929171,
        164.930319, 161.928775, 163.929197, 165.930290, 166.932045,
        167.932368, 169.935460, 168.934211, 167.933894, 169.934759,
        170.936322, 171.936378, 172.938207, 173.938858, 175.942568,
        174.940768, 175.942682, 173.940040, 175.941402, 176.943220,
        177.943698, 178.945815, 179.946549, 179.947466, 180.947996,
        179.946706, 181.948206, 182.950224, 183.950933, 185.954362,
        184.952956, 186.955751, 183.952491, 185.953838, 186.955748,
        187.955836, 188.958145, 189.958445, 191.961479, 190.960591,
        192.962924, 189.959930, 191.961035, 193.962664, 194.964774,
        195.964935, 197.967876, 196.966552, 195.965815, 197.966752,
        198.968262, 199.968309, 200.970285, 201.970626, 203.973476,
        202.972329, 204.974412, 203.973029, 205.974449, 206.975881,
        207.976636, 208.980383, 232.038050, 234.040946, 235.043923,
        238.050783,
        ])

    hydrogen_slice = slice(0, 2)
    helium_slice = slice(2, 4)
    metal_slice = slice(4, None)

    _check_default = True

    @abstractmethod
    def _abu_molfrac_raw(self, x):
        raise NotImplementedError()

    @abstractmethod
    def _abu_massfrac_raw(self, x):
        raise NotImplementedError()

    @abstractmethod
    def _abu_numfrac_raw(self, x):
        raise NotImplementedError()

    def _abu(self,
             x = 0.,
             molfrac = None,
             numfrac = None,
             normalize = True,
             ):
        """
        x is scale parameter ~ 'metallicity' approximator

        mass and mol fractions on input are normalized to
        sum(massfrac) == 1 if normalization is requested (default)

        number fractions are always normalized to 1 unless input is
        number fractions; then it is only normalized if requested (default)

        """

        if numfrac is None:
            numfrac = self.numfrac
        if molfrac is None:
            molfrac = self.molfrac

        try:
            abu = self._abu_molfrac_raw(x)

            if numfrac:
                if normalize:
                    abu /= np.sum(abu * self.masses)
                abu /= np.sum(abu)
            elif not molfrac:
                abu *= self.masses
                if normalize:
                    abu /= np.sum(abu)
            elif normalize:
                abu /= np.sum(abu * self.masses)
        except NotImplementedError:
            try:
                abu = self._abu_massfrac_raw(x)
                if numfrac:
                    abu /= self.masses
                if normalize or numfrac:
                    abu /= np.sum(abu)
                if molfrac and not numfrac:
                    abu /= self.masses
            except NotImplementedError:
                abu = self._abu_numfrac_raw(x)
                if numfrac:
                    if normalize:
                        abu /= np.sum(abu)
                else:
                    abu *= self.masses
                    abu /= np.sum(abu)
                    if molfrac:
                        abu /= self.masses
        return abu


    def _abu_metallicity(self, abu, Z_slice = None):
        if Z_slice is None:
            Z_slice = self.metal_slice
        return np.sum(abu[Z_slice])

    def _abu_hydrogen(self, abu):
        return np.sum(abu[self.hydrogen_slice])

    def _abu_z(self, abu):
        z = self._abu_metallicity(abu, Z_slice = self.Z_slice) / self.solar_metal
        if self.href:
            h_abu = self._abu_metallicity(abu, Z_slice = self.hydrogen_slice)
            z *= self.h_sun / h_abu
        return z

    def _z_xi(self, x):
        return self._abu_z(self._abu(x))

    def _find_x(self, z0, **kw):
        if z0 == 0:
            return z0
        # start assume xi == z as a starting guess
        x = z0

        # setup
        z = self._z_xi(x)
        r = np.array([[z, x]]*2, dtype = np.float)

        while abs(1 - z/z0) > 1.e-12:
            if r[0,0] > z0:
                x = r[0,1] * 0.5
                z = self._z_xi(x)
                r[1,:] = r[0,:]
                r[0,:] = [z, x]
            elif r[1,0] < z0:
                x = r[1,1] * 2
                z = self._z_xi(x)
                r[0,:] = r[1,:]
                r[1,:] = [z, x]
            else:
                # simple for now
                # x = 0.5 * ( r[0,1] + r[1,1] )

                # better - linear interpolation
                f = (z0 - r[0,0]) / (r[1,0] - r[0,0])
                x = (1-f) * r[0,1] + f * r[1,1]

                z = self._z_xi(x)
                if z > z0:
                    r[1,:] = [z, x]
                else:
                    r[0,:] = [z, x]
        return x

    def __init__(self,
                 z = 1,
                 silent = False,
                 check = None,
                 molfrac = False,
                 numfrac = False,
                 massfrac = None,
                 metal = None,
                 href = False,
                 ):
        """
        Generate abundances set using z = Z/Z_sun, the scale relative
        to solar, if z is positive; absolute metallicity Z = -z if z
        is negative.

        molfrac = metallicity my mol fraction

        numfrac = metallicity my number fraction

        massfrac = metallicity my mass fraction [default]

        metal = normalize to given metal, e.g., 'Fe'

        href = use fraction relative to hydrogen, e.g., to compute [Fe/H]

        """
        super().__init__()

        self.setup_logger(silent = silent)

        if check is None:
            check = self._check_default
        if check:
            numiso = len(self.ions)
            # assert self.ions == IonList(SolAbu('Lo09').ions())
            assert self.masses.shape[0] == numiso

        if massfrac is True:
            assert molfrac == numfrac == False
        else:
            massfrac = not (molfrac or numfrac)
        assert np.count_nonzero([massfrac, molfrac, numfrac]) == 1

        if metal is not None:
            ion_ref = ion(metal)
            Z_slice = np.where(ufunc_Z(np.array(self.ions)) == ion_ref.Z)[0]
        else:
            Z_slice = None
        self.Z_slice = Z_slice

        self.numfrac = numfrac
        self.molfrac = molfrac
        self.href = href

        self.solar_abu  = self._abu(1.)
        self.solar_metal  = self._abu_metallicity(self.solar_abu, Z_slice = self.Z_slice)
        if href:
            self.h_sun = self._abu_metallicity(self.solar_abu, Z_slice = self.hydrogen_slice)

        # negative values are interpreted as absolute values of z
        # rather than scale factor
        if z < 0:
            z = -z / self.solar_metal

        x = self._find_x(z)

        abu = self._abu(x, molfrac = False, numfrac = False)

        self.abu = abu
        self.iso = np.array(self.ions)

        self.normalize()
        self.is_sorted = True
        self.sort()

        self.close_logger(timing='Abundance set generated in')
