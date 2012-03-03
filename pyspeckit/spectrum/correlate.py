import numpy as np
import arithmetic
import units
import classes
import pyfits
import headers

def correlate(spectrum1, spectrum2, range=None, units=None):
    """
    Cross-correlate spectrum1 with spectrum2

    """

    if range is not None:
        spectrum1 = spectrum1.slice(*range, units=units)
        spectrum2 = spectrum2.slice(*range, units=units)

    if not all(spectrum1.xarr1 == spectrum2.xarr2):
        spectrum2 = arithmetic.interp(spectrum2, spectrum1)

    data1 = spectrum1.data
    data2 = spectrum2.data

    xcorr = np.correlate(data1, data2, mode='same')
    
    xarr = spectrum1.xarr
    xrange = xarr.max()-xarr.min()
    xmin = -xrange/2.
    xmax =  xrange/2.
    offset_values = np.linspace(xmin, xmax, len(xarr))

    offset_xarr = units.SpectroscopicAxis(offset_values, unit=xarr.units) 

    header = headers.intersection(spectrum1.header, spectrum2.header)
    header['CRPIX1'] = 1
    header['CRVAL1'] = xmin
    header['CDELT1'] = offset_xarr.cdelt()

    return classes.XCorrSpectrum(xarr=offset_xarr, data=xcorr, header=header)

