import numpy as np
import interpolation
import units as units_module
import classes
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
import headers

def correlate(spectrum1, spectrum2, range=None, units=None, errorweight=False):
    """
    Cross-correlate spectrum1 with spectrum2

    """

    if range is not None:
        spectrum1 = spectrum1.slice(*range, units=units)
        spectrum2 = spectrum2.slice(*range, units=units)

    if not (spectrum1.xarr.shape == spectrum2.xarr.shape) or not all(spectrum1.xarr == spectrum2.xarr):
        spectrum2 = interpolation.interp(spectrum2, spectrum1)

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
    error = (np.correlate( (spectrum1.error)**2 , np.ones(xcorr.shape), mode='same') + 
             np.correlate( (spectrum2.error)**2 , np.ones(xcorr.shape), mode='same'))**0.5
    
    xarr = spectrum1.xarr
    xrange = xarr.max()-xarr.min()
    xmin = -xrange/2.
    xmax =  xrange/2.
    offset_values = np.linspace(xmin, xmax, len(xarr))

    offset_xarr = units_module.SpectroscopicAxis(offset_values, unit=xarr.units) 

    header = headers.intersection(spectrum1.header, spectrum2.header)
    header.update('CRPIX1',1)
    try:
        header.update('CRVAL1',xmin)
    except ValueError:
        header.update('CRVAL1',xmin.tolist())
    header.update('CDELT1',offset_xarr.cdelt())

    return classes.XCorrSpectrum(xarr=offset_xarr, data=xcorr, header=header, error=error)

