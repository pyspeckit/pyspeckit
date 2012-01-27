from .. import units
import numpy as np
import numpy.ma as ma
import read_class

readers = {}
suffix_types = {}

def make_axis(xarr,hdr,specname=None, wcstype='', specaxis="1", verbose=True):
    """
    Parse parameters from a .fits header into required SpectroscopicAxis
    parameters
    """

    #DEBUG if wcstype is not '': print "Loading file with WCSTYPE %s" % wcstype

    xunits = hdr.get('CUNIT%s%s' % (specaxis,wcstype))
    if hdr.get('ORIGIN') == 'CLASS-Grenoble' and xunits is None:
        # CLASS default
        xunits = 'Hz'
        
    # SDSS doesn't use FITS standard! Argh.    
    if hdr.get('TELESCOP') == 'SDSS 2.5-M':   
        xunits = 'angstroms' 

    # IRAF also doesn't use the same standard
    if xunits is None:
        if hdr.get('WAT1_001') is not None:
            pairs = hdr.get('WAT1_001').split()
            pdict = dict( [s.split("=") for s in pairs] )
            if 'units' in pdict:
                xunits = pdict['units']

    if hdr.get('REFFREQ'+wcstype):
        refX = hdr.get('REFFREQ'+wcstype)
    elif hdr.get('RESTFREQ'+wcstype):
        refX = hdr.get('RESTFREQ'+wcstype)
    elif hdr.get('RESTFRQ'+wcstype):
        refX = hdr.get('RESTFRQ'+wcstype)
    else:
        if verbose: print "No reference frequency found.  Velocity transformations will not be possible"
        refX = None

    if hdr.get('CTYPE%s%s' % (specaxis,wcstype)):
        xtype = hdr.get('CTYPE%s%s' % (specaxis,wcstype))
    else:
        xtype = 'VLSR'

    XAxis = units.SpectroscopicAxis(xarr,xunits,xtype=xtype,refX=refX)

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

import fits_reader
open_1d_fits = check_reader(fits_reader.open_1d_fits)
open_1d_pyfits = check_reader(fits_reader.open_1d_pyfits)
import tspec_reader
tspec_reader = check_reader(tspec_reader.tspec_reader)
import txt_reader
open_1d_txt = check_reader(txt_reader.open_1d_txt)
import hdf5_reader
open_hdf5 = check_reader(hdf5_reader.open_hdf5)
