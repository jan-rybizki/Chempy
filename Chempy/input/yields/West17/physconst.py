import math

dnucs133hfs = 9192631770 # ?: (genauer ??(133Cs)hfs) Casiumfrequenz, exakt 9192631770 Hz, die Frequenz der Strahlung beim ?bergang zwischen den beiden Hyperfeinstrukturniveaus des Grundzustandes von Atomen des Nuklids 133Cs.
CLIGHT = 29979245800 # : Lichtgeschwindigkeit im Vakuum, exakt 299792458 Meter/Sekunde
PLANCK = 6.62606957e-27 # h: Planck'sches Wirkungsquantum, exakt 6,62606957 * 10-34 Joule?Sekunde
ECHARGE = 1.602176565e-19 * CLIGHT / 10 # e: Elementarladung, exakt 1,602176565 * 10-19 Coulomb,
KB = 1.3806488e-16 # k : Boltzmann-Konstante, exakt 1,3806488 * 10-23 Joule/Kelvin,
NA = 6.02214129e23 # NA: Avogadro-Konstante, exakt 6,02214129 * 10+23 1/mol
KCD = 683 # Kcd: Lichtausbeute, exakt 683 Lumen/Watt bei monochromatischer Strahlung mit 540 * 1012 Hertz





SEC     = 31556926.
PI      = 4 * math.atan(1.)
XMSUN   = 1.9891e33
XLSUN   = 3.84e33
RSUN    = 6.98e10
GRAV    = 6.67259e-8
SB      = 5.67051e-5
# KB      = 1.380658e-16
# PLANCK  = 6.6260755e-27
HBAR    = 1.05457266e-27
# CLIGHT  = 29979245888e0
PC      = 3.08567756706e18
AU      = 1.4959787e13
ECHARGE = 4.8032067993e-10
ME      = 9.1093897e-28
EV      = 1.60217733e-12
MEV     = EV * 1.e6
RK      = 83145112.1195e0
ARAD    = 7.56591414985e-15
# NA      = 6.0221367e23
AMU     = 1.6605402e-24
DAY     = 86400.0e0
BARN    = 1.e-24
SIGT    = 8 * PI / 3 * (ECHARGE**2 / (ME * CLIGHT**2))**2

class Kepler:
    """
    Physics constants as used with KEPLER
    """
    solmass  = 1.9892e+33
    gee      = 6.670e-8
    n0       = 6.02254e+23
    sigt     = 6.65205e-25
    k        = 1.38054e-16
    a        = 7.5648e-15
    me       = 9.10908e-28
    h        = 6.62559e-27
    c        = 2.997925e+10
    pie      = 3.1415926536e+0
    solmassi = 1 / solmass
    solrad   = 6.9599e+10
    solradi  = 1 / solrad
    penmex   = 0.782333098e0 # (m_n-m_p-m_e)*c**2 [MeV]
    year     = 3.1558e+7
    pie43    = pie * 4 / 3

    rk       = k * n0
    sb       = a * c / 4
