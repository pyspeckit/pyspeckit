import pyfits
from .. import units
import numpy.ma as ma
import numpy as np
from . import make_axis
import operator

def open_1d_fits(filename,**kwargs):
    """
    Grabs all the relevant pieces of a simple FITS-compliant 1d spectrum

    Inputs:
        wcstype - the suffix on the WCS type to get to
            velocity/frequency/whatever
        specnum - Which # spectrum, along the y-axis, is 
            the data?
        errspecnum - Which # spectrum, along the y-axis,
            is the error spectrum?

    """

    f = pyfits.open(filename)

    return open_1d_pyfits(f[0],**kwargs)


def open_1d_pyfits(pyfits_hdu,specnum=0,wcstype='',specaxis="1",errspecnum=None, autofix=True,
        scale_keyword=None, scale_action=operator.div, **kwargs):
    """
    This is open_1d_fits but for a pyfits_hdu so you don't necessarily have to
    open a fits file
    """

    hdr = pyfits_hdu._header
    if autofix: 
        for card in hdr.ascardlist():
            try:
                card.verify('silentfix')
            except pyfits.VerifyError:
                hdr.__delitem__(card.key)
    data = pyfits_hdu.data

    # search for the correct axis (may be 1 or 3, unlikely to be 2 or others)
    # 1 = 1D spectrum
    # 3 = "3D" spectrum with a single x,y point (e.g., JCMT smurf/makecube)
    if hdr.get('NAXIS') > 1:
        for ii in xrange(1,hdr.get('NAXIS')+1):
            ctype = hdr.get('CTYPE%i'%ii)
            if ctype in units.xtype_dict:
                specaxis="%i" % ii

    if hdr.get('NAXIS') == 2:
        if isinstance(specnum,list):
            # allow averaging of multiple spectra (this should be modified
            # - each spectrum should be a Spectrum instance)
            spec = ma.array(data[specnum,:]).mean(axis=0)
        elif isinstance(specnum,int):
            spec = ma.array(data[specnum,:]).squeeze()
        else:
            raise TypeError("Specnum is of wrong type (not a list of integers or an integer).  Type: %s" %
                    str(type(specnum)))
        if errspecnum is not None:
            errspec = ma.array(data[errspecnum,:]).squeeze()
        else:
            errspec = spec*0 # set error spectrum to zero if it's not in the data

    elif hdr.get('NAXIS') > 2:
        for ii in xrange(2,hdr.get('NAXIS')):
            # only fail if extra axes have more than one row
            if hdr.get('NAXIS%i' % ii) > 1:
                raise ValueError("Too many axes for open_1d_fits")
        spec = ma.array(data).squeeze()
        if errspecnum is None: errspec = spec*0 # set error spectrum to zero if it's not in the data
    else:
        spec = ma.array(data).squeeze()
        if errspecnum is None: errspec = spec*0 # set error spectrum to zero if it's not in the data

    if hdr.get(scale_keyword):
        print "Found SCALE keyword %s.  Using %s to scale it" % (scale_keyword,scale_action)
        scaleval = hdr.get(scale_keyword)
        spec = scale_action(spec,scaleval)
        errspec = scale_action(errspec,scaleval)

    xarr = None
    if hdr.get('ORIGIN') == 'CLASS-Grenoble':
        # Use the CLASS FITS definition (which is non-standard)
        # http://iram.fr/IRAMFR/GILDAS/doc/html/class-html/node84.html
        # F(n) = RESTFREQ + CRVALi + ( n - CRPIXi ) * CDELTi
        print "Loading a CLASS .fits spectrum"
        dv = -1*hdr.get('CDELT1')
        v0 = hdr.get('RESTFREQ') + hdr.get('CRVAL1')
        p3 = hdr.get('CRPIX1')
    elif hdr.get('CD%s_%s%s' % (specaxis,specaxis,wcstype)):
        dv,v0,p3 = hdr['CD%s_%s%s' % (specaxis,specaxis,wcstype)],hdr['CRVAL%s%s' % (specaxis,wcstype)],hdr['CRPIX%s%s' % (specaxis,wcstype)]
        hdr.update('CDELT%s' % specaxis,dv)
    elif hdr.get('CDELT%s%s' % (specaxis,wcstype)):
        dv,v0,p3 = hdr['CDELT%s%s' % (specaxis,wcstype)],hdr['CRVAL%s%s' % (specaxis,wcstype)],hdr['CRPIX%s%s' % (specaxis,wcstype)]
    elif len(data.shape) > 1:
        # try assuming first axis is X axis
        if hdr.get('CUNIT%s%s' % (specaxis,wcstype)):
            xarr = data[0,:]
            spec = data[1,:]
            if data.shape[0] > 2:
                errspec = data[2,:]
        else:
            raise TypeError("Don't know what type of FITS file you've input; its header is not FITS compliant and it doesn't look like it was written by pyspeckit.")

    # Deal with logarithmic wavelength binning if necessary
    if xarr is None:
        if hdr.get('WFITTYPE') == 'LOG-LINEAR':
            xconv = lambda v: 10**((v-p3+1)*dv+v0)
            xarr = xconv(np.arange(len(spec)))
        else:
            xconv = lambda v: ((v-p3+1)*dv+v0)
            xarr = xconv(np.arange(len(spec)))
    
    restfreq = hdr.get('RESTFREQ')
    if restfreq is None: restfrq= hdr.get('RESTFRQ')

    XAxis = make_axis(xarr,hdr,wcstype=wcstype,specaxis=specaxis,**kwargs)

    return spec,errspec,XAxis,hdr

