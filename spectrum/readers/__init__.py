import spectrum.units as units
import numpy as np
import numpy.ma as ma

readers = {}
suffix_types = {}

def make_axis(xarr,hdr,specname=None, wcstype=''):
    """
    Parse parameters from a .fits header into required SpectroscopicAxis
    parameters
    """

    xunits = hdr.get('CUNIT1'+wcstype)
    if hdr.get('ORIGIN') == 'CLASS-Grenoble' and xunits is None:
        # CLASS default
        xunits = 'Hz'

    if hdr.get('REFFREQ'+wcstype):
        reffreq = hdr.get('REFFREQ'+wcstype)
    elif hdr.get('RESTFREQ'+wcstype):
        reffreq = hdr.get('RESTFREQ'+wcstype)
    elif hdr.get('RESTFRQ'+wcstype):
        reffreq = hdr.get('RESTFRQ'+wcstype)
    else:
        reffreq = None

    if hdr.get('CTYPE1'+wcstype):
        xtype = hdr.get('CTYPE1'+wcstype)
    else:
        xtype = 'VLSR'

    XAxis = units.SpectroscopicAxis(xarr,xunits,xtype=xtype,reffreq=reffreq)

    return XAxis

from fits_reader import open_1d_fits
from tspec_reader import tspec_reader
from txt_reader import open_1d_txt
