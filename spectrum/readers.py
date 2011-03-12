import numpy as np
import numpy.ma as ma

def open_1d_txt(filename):
    import atpy
    T = atpy.Table(filename, type='ascii',
            Reader=atpy.asciitables.asciitable.CommentedHeader, masked=True)

    
    xarr = T.data[T.data.dtype.names[0]]
    data = T.data[T.data.dtype.names[1]]
    if len(T.columns) > 2:
        error = T.data[T.data.dtype.names[2]]
    else:
        # assume uniform, nonzero error
        error = data*0 + 1.0

    return xarr,data,error,T


def open_1d_fits(filename,specnum=0,wcstype='',errspecnum=None,**kwargs):
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
    import pyfits
    f = pyfits.open(filename)
    hdr = f[0].header
    spec = ma.array(f[0].data).squeeze()
    errspec  = None
    if hdr.get('NAXIS') == 2:
        if errspecnum is not None:
            errspec = spec[errspecnum,:]
        if isinstance(specnum,list):
            # allow averaging of multiple spectra (this should be modified
            # - each spectrum should be a Spectrum instance)
            spec = spec[specnum,:].mean(axis=0)
        elif isinstance(specnum,int):
            spec = spec[specnum,:]
        else:
            raise TypeError("Specnum is of wrong type (not a list of integers or an integer).  Type: %s" %
                    str(type(specnum)))
    elif hdr.get('NAXIS') > 2:
        for ii in xrange(2,hdr.get('NAXIS')):
            # only fail if extra axes have more than one row
            if hdr.get('NAXIS%i' % ii) > 1:
                raise ValueError("Too many axes for open_1d_fits instead")
    if hdr.get('ORIGIN') == 'CLASS-Grenoble':
        # Use the CLASS FITS definition (which is non-standard)
        # http://iram.fr/IRAMFR/GILDAS/doc/html/class-html/node84.html
        # F(n) = RESTFREQ + CRVALi + ( n - CRPIXi ) * CDELTi
        print "Loading a CLASS .fits spectrum"
        dv = -1*hdr.get('CDELT1')
        v0 = hdr.get('RESTFREQ') + hdr.get('CRVAL1')
        p3 = hdr.get('CRPIX1')
    elif hdr.get('CD1_1'+wcstype):
        dv,v0,p3 = hdr['CD1_1'+wcstype],hdr['CRVAL1'+wcstype],hdr['CRPIX1'+wcstype]
    else:
        dv,v0,p3 = hdr['CDELT1'+wcstype],hdr['CRVAL1'+wcstype],hdr['CRPIX1'+wcstype]

    # Deal with logarithmic wavelength binning if necessary
    if hdr.get('WFITTYPE') == 'LOG-LINEAR':
        xconv = lambda v: 10**((v-p3+1)*dv+v0)
        xarr = xconv(np.arange(len(spec)))
    else:
        xconv = lambda v: ((v-p3+1)*dv+v0)
        xarr = xconv(np.arange(len(spec)))

    XAxis = make_axis(xarr,hdr,**kwargs)

    return spec,errspec,XAxis,hdr

def make_axis(xarr,hdr,specname=None, wcstype=''):
    """
    Parse parameters from a .fits header into required SpectroscopicAxis
    parameters
    """
    import units

    xunits = hdr.get('CUNIT1'+wcstype)
    if hdr.get('ORIGIN') == 'CLASS-Grenoble' and xunits is None:
        # CLASS default
        xunits = 'Hz'

    if hdr.get('REFFREQ'+wcstype):
        reffreq = hdr.get('REFFREQ'+wcstype)
    elif hdr.get('RESTFREQ'+wcstype):
        reffreq = hdr.get('RESTFREQ'+wcstype)
    else:
        reffreq = None

    if hdr.get('CTYPE1'+wcstype):
        xtype = hdr.get('CTYPE1'+wcstype)
    else:
        xtype = 'VLSR'

    XAxis = units.SpectroscopicAxis(xarr,xunits,xtype=xtype,reffreq=reffreq)

    return XAxis
