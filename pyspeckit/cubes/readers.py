"""
Author: Adam Ginsburg
Created: 3/17/2011
"""
from __future__ import print_function
import numpy as np
import numpy.ma as ma
from .. import spectrum
import operator
from astropy import log

def open_3d_fits(filename,wcstype='',average_extra=False, specaxis=3,
                 scale_keyword=None, scale_action=operator.div, **kwargs):
    """
    Grabs all the relevant pieces of a simple FITS-compliant 3d data cube

    Parameters
    ----------
    wcstype : str
        the suffix on the WCS type to get to
        velocity/frequency/whatever
    specaxis : int 
        Which axis containts the spectrum?  Default 3
    scale_keyword : str
        A string to use to scale the data using the action
        ``scale_action``

    """
    try:
        import astropy.io.fits as pyfits
    except ImportError:
        import pyfits
    if isinstance(filename, pyfits.hdu.image.PrimaryHDU):
        f=[filename]
    elif isinstance(filename, pyfits.hdu.hdulist.HDUList):
        f=filename
    else:
        f = pyfits.open(filename,ignore_missing_end=True)
    hdr = f[0].header
    cube = ma.array(f[0].data.squeeze())  # remove extra dimensions such as polarization
    cube.mask = np.isnan(cube) + np.isinf(cube)
    if hdr.get('NAXIS') > 3 and average_extra is False:
        print("Extra dimensions in data cube will be ignored (only the first element of that dimension will be kept).")
    while len(cube.shape) > 3:
        if average_extra:
            cube = cube.mean(axis=0)
        else:
            cube = cube[0,...].squeeze()

    # Does CLASS even make cubes?  Maybe....
    if hdr.get('ORIGIN') == 'CLASS-Grenoble':
        # Use the CLASS FITS definition (which is non-standard)
        # http://iram.fr/IRAMFR/GILDAS/doc/html/class-html/node84.html
        # F(n) = RESTFREQ + CRVALi + ( n - CRPIXi ) * CDELTi
        print("Loading a CLASS .fits spectrum")
        dv = -1*hdr.get('CDELT3')
        v0 = hdr.get('RESTFREQ') + hdr.get('CRVAL3')
        p3 = hdr.get('CRPIX3')
    elif hdr.get(('CD%i_%i' % (specaxis,specaxis))+wcstype):
        dv,v0,p3 = hdr[('CD%i_%i' % (specaxis,specaxis))+wcstype],hdr[('CRVAL%i' % specaxis)+wcstype],hdr[('CRPIX%i' % specaxis)+wcstype]
    else:
        dv,v0,p3 = hdr[('CDELT%i' % specaxis)+wcstype],hdr[('CRVAL%i' % specaxis)+wcstype],hdr[('CRPIX%i' % specaxis)+wcstype]

    if scale_keyword is not None:
        if scale_keyword in hdr:
            scaleval = hdr[scale_keyword]
            log.info("Found SCALE keyword %s.  Using %s to scale it" % (scale_keyword,scale_action))
            cube = scale_action(cube,scaleval)
        else:
            raise KeyError('{0} not found in header for file {1}'.format(scale_keyword, filename))

    xconv = lambda v: ((v-p3+1)*dv+v0)
    xarr = xconv(np.arange(cube.shape[3-specaxis]))

    XAxis = spectrum.readers.make_axis(xarr,hdr,wcstype=wcstype, specaxis=specaxis, **kwargs)

    return cube,XAxis,hdr,f
