#! /bin/env python3

"""
Module for human-readable byte strings.
"""

import sys

_units = ('','k','M','G','T','P','E','Z','Y')

def byte2human(size,
               strip = True,
               SI = False,
               short = False,
               length = 3):
    """
    Return byte string in human-readible format.
    """
    assert 2 < length < 6, 'Length out of Range'
    su = 'B'
    if SI:
        prefix = ''
        div = 1000
    else:
        prefix = 'i'
        div = 1024
    div_lim = 1000 * (1 - 0.5 * 10**(-length) - 2.e-15)
    if short:
        prefix = ''
        su = ''

    asize = abs(size)
    osize = round(asize)
    assert abs((osize + 1.e-16)/(asize + 1.e-16) - 1) < 1.e-14
    xsize = osize

    i = 0
    while xsize > div_lim:
        xsize /= div
        i += 1
    if length == 3:
        rnd_lim = 10 * (1 - 0.5 * 10**(-length + 1) - 2.e-15)
        if (xsize >= rnd_lim) or (i == 0):
            sv = '{:3d}'.format(int(round(xsize)))
        else:
            sv = '{:3.1f}'.format(xsize)
    elif length == 4:
        rnd_lim1 = 100 * (1 - 0.5 * 10**(-length + 1) - 2.e-15)
        rnd_lim2 = 10 * (1 - 0.5 * 10**(-length + 1) - 2.e-15)
        if (xsize >= rnd_lim1) or (i == 0):
            sv = '{:4d}'.format(int(round(xsize)))
        elif (xsize >= rnd_lim2):
            sv = '{:4.1f}'.format(xsize)
        else:
            sv = '{:4.2f}'.format(xsize)
    elif length == 5:
        rnd_lim1 = 1000 * (1 - 0.5 * 10**(-length + 1) - 2.e-15)
        rnd_lim2 = 100 * (1 - 0.5 * 10**(-length + 1) - 2.e-15)
        rnd_lim3 = 10 * (1 - 0.5 * 10**(-length + 1) - 2.e-15)
        if i > 0 and 999 < round(xsize*div) < div:
            xsize *= div
            i -= 1
            sv = '{:4d}'.format(int(round(xsize)))
            sv = sv[0] + ',' + sv[1:4]
        elif (xsize >= rnd_lim1) or (i == 0):
            sv = '{:5d}'.format(int(round(xsize)))
        elif (xsize >= rnd_lim2):
            sv = '{:5.1f}'.format(xsize)
        elif (xsize >= rnd_lim3):
            sv = '{:5.2f}'.format(xsize)
        else:
            sv = '{:5.3f}'.format(xsize)
    else:
        raise Exception('Length out of Range')

    if i >= len(_units):
        sv = '*'*length
    unit = _units[i]
    if i >= 1:
        unit += prefix
    unit = unit + su

    if round(size) < 0:
        sv = '-' + sv.strip()
        sv = ' '*(length - len(sv)) + sv
    s = sv + ' ' + unit
    if strip:
        s = s.strip()

    return s


if __name__ == '__main__':
    argv = sys.argv
    if len(argv) == 2:
        try:
            print(byte2human(round(float(argv[1]))))
        except:
            print('***')
