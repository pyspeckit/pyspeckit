from __future__ import print_function
import numpy as np
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
from . import interpolation
from . import units as units_module
from . import classes
from . import headers

def correlate(spectrum1, spectrum2, srange=None, unit=None, errorweight=False):
    """
    Cross-correlate spectrum1 with spectrum2

    """
    def my_all(barr):
        """ why does all on a boolean array not work???
        
            TypeError: 'bool' object is not iterable

        """
        print("my_all:",type(barr))     # my_all: <class 'bool'>
        for b in barr:
            if not b: return False
        return True

    if srange is not None:
        spectrum1 = spectrum1.slice(*srange, unit=unit)
        spectrum2 = spectrum2.slice(*srange, unit=unit)

    if True:
        # the all() does not work anymore 
        # PJT -> TypeError: 'bool' object is not iterable
        print(spectrum1.xarr.shape)
        print(spectrum2.xarr.shape)
        print(spectrum1.xarr)
        print(spectrum2.xarr)
        is_1 = spectrum1.xarr.shape == spectrum2.xarr.shape
        print(type(is_1),is_1)
        is_2 = spectrum1.xarr == spectrum2.xarr
        print(type(is_2),is_2)
        is_3 = my_all(is_2)
        #if not (spectrum1.xarr.shape == spectrum2.xarr.shape) or not all(spectrum1.xarr == spectrum2.xarr):
        if not is_1 or not is_3:
            spectrum2 = interpolation.interp(spectrum2, spectrum1)
        else:
            print("PJT: no interpolation needed")
    else:
        print("PJT: interpolation disabled")


    data1 = spectrum1.data
    data2 = spectrum2.data

    xcorr = np.correlate(data1, data2, mode='same')

    # very simple propagation of error
    # each element is multiplied, multiplicative error is given such that (sigma_xy/xy)**2 = (sigma_x/x)**2 + (sigma_y/y)**2
    # error = (np.correlate( (spectrum1.error/spectrum1.data)**2 , np.ones(xcorr.shape), mode='same') +
    #          np.correlate( (spectrum2.error/spectrum2.data)**2 , np.ones(xcorr.shape), mode='same'))**0.5 * xcorr
    # That approach sucks - what if data == 0?
    #
    # this might be more correct: http://arxiv.org/pdf/1006.4069v1.pdf eqn 4
    # but it doesn't quite fit my naive expectations so:
    error = (np.correlate((spectrum1.error)**2, np.ones(xcorr.shape), mode='same') +
             np.correlate((spectrum2.error)**2, np.ones(xcorr.shape), mode='same'))**0.5

    xarr = spectrum1.xarr
    x_range = xarr.max()-xarr.min()
    xmin = -x_range/2.
    xmax =  x_range/2.
    offset_values = np.linspace(xmin, xmax, len(xarr))

    offset_xarr = units_module.SpectroscopicAxis(offset_values, unit=xarr.unit)

    header = headers.intersection(spectrum1.header, spectrum2.header)
    header['CRPIX1'] = 1
    try:
        header['CRVAL1'] = xmin
    except ValueError:
        try:
            header['CRVAL1'] = xmin.tolist()
        except NotImplementedError:
            header['CRVAL1'] = xmin.value
    try:
        header['CDELT1'] = offset_xarr.cdelt()
    except ValueError:
        header['CDELT1'] = offset_xarr.cdelt().value

    return classes.XCorrSpectrum(xarr=offset_xarr, data=xcorr, header=header, error=error)
