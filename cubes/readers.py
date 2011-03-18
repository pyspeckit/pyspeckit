"""
Author: Adam Ginsburg
Created: 3/17/2011
"""
import numpy as np
import numpy.ma as ma
import spectrum

def open_3d_fits(filename,wcstype='',average_extra=False, **kwargs):
    """
    Grabs all the relevant pieces of a simple FITS-compliant 3d data cube

    Inputs:
        wcstype - the suffix on the WCS type to get to
            velocity/frequency/whatever

    """
    import pyfits
    f = pyfits.open(filename)
    hdr = f[0].header
    cube = ma.array(f[0].data).squeeze()  # remove extra dimensions such as polarization
    cube.mask = np.isnan(cube) + np.isinf(cube)
    if hdr.get('NAXIS') > 3 and average_extra is False:
        print "Extra dimensions in data cube will be ignored (only the first element of that dimension will be kept)."
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
        print "Loading a CLASS .fits spectrum"
        dv = -1*hdr.get('CDELT3')
        v0 = hdr.get('RESTFREQ') + hdr.get('CRVAL3')
        p3 = hdr.get('CRPIX3')
    elif hdr.get('CD3_3'+wcstype):
        dv,v0,p3 = hdr['CD3_3'+wcstype],hdr['CRVAL3'+wcstype],hdr['CRPIX3'+wcstype]
    else:
        dv,v0,p3 = hdr['CDELT3'+wcstype],hdr['CRVAL3'+wcstype],hdr['CRPIX3'+wcstype]

    xconv = lambda v: ((v-p3+1)*dv+v0)
    xarr = xconv(np.arange(cube.shape[0]))

    XAxis = spectrum.readers.make_axis(xarr,hdr,wcstype=wcstype,**kwargs)

    return cube,XAxis,hdr,f
