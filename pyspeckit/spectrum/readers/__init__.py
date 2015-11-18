from __future__ import print_function
import numpy as np
from astropy import units as u

from .. import units
from ...specwarnings import warn

readers = {}
suffix_types = {}

def _parse_velocity_convention(vc):
    if vc in (u.doppler_radio, 'radio', 'RADIO', 'VRAD', 'F', 'FREQ'):
        return 'radio'
    elif vc in (u.doppler_optical, 'optical', 'OPTICAL', 'VOPT', 'W', 'WAVE'):
        return 'optical'
    elif vc in (u.doppler_relativistic, 'relativistic', 'RELATIVE', 'VREL',
                'speed', 'V', 'VELO'):
        return 'relativistic'

def make_axis(xarr, hdr, specname=None, wcstype='', specaxis="1", verbose=True,
              **kwargs):
    """
    Parse parameters from a .fits header into required SpectroscopicAxis
    parameters
    """

    # DEBUG if wcstype is not '': print "Loading file with WCSTYPE %s" % wcstype

    xunits = hdr.get('CUNIT%s%s' % (specaxis,wcstype))
    if hdr.get('ORIGIN') == 'CLASS-Grenoble' and xunits is None:
        # CLASS default
        xunits = 'Hz'
        
    # SDSS doesn't use FITS standard! Argh.
    if hdr.get('TELESCOP') == 'SDSS 2.5-M':
        xunits = 'angstrom'

    # IRAF also doesn't use the same standard
    if xunits is None:
        if hdr.get('WAT1_001') is not None:
            pairs = hdr.get('WAT1_001').split()
            pdict = dict([s.split("=") for s in pairs])
            if 'units' in pdict:
                xunits = pdict['units']

    if hdr.get('REFFREQ'+wcstype):
        refX = hdr.get('REFFREQ'+wcstype)*u.Hz
    elif hdr.get('RESTFREQ'+wcstype):
        refX = hdr.get('RESTFREQ'+wcstype)*u.Hz
    elif hdr.get('RESTFRQ'+wcstype):
        refX = hdr.get('RESTFRQ'+wcstype)*u.Hz
    else:
        if verbose:
            warn("Warning: No reference frequency found."
                 "  Velocity transformations will not be "
                 "possible unless you set a reference frequency/wavelength")
        refX = None

    if hdr.get('VELDEF'):
        convention, frame = units.parse_veldef(hdr['VELDEF'])
        # vframe = hdr.get('VFRAME') if hdr.get('VFRAME') is not None else 0.0
    else:
        convention, frame = _parse_velocity_convention(hdr.get('CTYPE%s%s' % (specaxis,wcstype))), None
        # vframe = 0.0

    XAxis = units.SpectroscopicAxis(xarr, xunits, refX=refX,
                                    velocity_convention=convention, **kwargs)

    return XAxis

class ReaderError(Exception):
    pass

def check_reader(func):
    def reader(*args,**kwargs):
        returns = func(*args,**kwargs)
        if len(returns) != 4:
            raise ReaderError("Error: reader returns %i parameters instead of 4." % len(returns))
        else:
            data,error,xarr,header = returns
            if data.shape != error.shape:
                raise ValueError("Data and error spectra shapes do not match.")
            if data.shape != xarr.shape:
                raise ValueError("Data and X-axis shapes do not match.")
        return returns
    return reader

from . import fits_reader
open_1d_fits = check_reader(fits_reader.open_1d_fits)
open_1d_pyfits = check_reader(fits_reader.open_1d_pyfits)
from . import tspec_reader
tspec_reader = check_reader(tspec_reader.tspec_reader)
from . import txt_reader
open_1d_txt = check_reader(txt_reader.open_1d_txt)
from . import hdf5_reader
open_hdf5 = check_reader(hdf5_reader.open_hdf5)
from .gbt import GBTSession
from .sdss_reader import read_sdss
from .galex import read_galex
