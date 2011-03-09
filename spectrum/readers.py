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
    spec = ma.array(f[0].data)
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
        raise ValueError("Too many axes for open_1d - use open_3d instead")
    if hdr.get('CD1_1'+wcstype):
        dv,v0,p3 = hdr['CD1_1'+wcstype],hdr['CRVAL1'+wcstype],hdr['CRPIX1'+wcstype]
    else:
        dv,v0,p3 = hdr['CDELT1'+wcstype],hdr['CRVAL1'+wcstype],hdr['CRPIX1'+wcstype]

    # Deal with logarithmic wavelength binning if necessary
    try: 
        if hdr['WFITTYPE'] == 'LOG-LINEAR':
            xconv = lambda v: 10**((v-p3+1)*dv+v0)
            xarr = xconv(np.arange(len(spec)))
    except KeyError:
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

    if hdr.get('REFFREQ'+wcstype):
        reffreq = hdr.get('REFFREQ'+wcstype)
    else:
        reffreq = None

    if hdr.get('CTYPE1'+wcstype):
        xtype = hdr.get('CTYPE1'+wcstype)
    else:
        xtype = 'VLSR'

    XAxis = units.SpectroscopicAxis(xarr,xunits,xtype=xtype,reffreq=reffreq)

    return XAxis
