import units
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

    if 'xunits' in T.keywords:
        xunits = T.keywords['xunits']
    else:
        xunits = 'unknown'

    XAxis = units.SpectroscopicAxis(xarr,xunits)

    return XAxis,data,error,T


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
    if errspecnum is None: errspec = spec*0 # set error spectrum to zero if it's not in the data
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
                raise ValueError("Too many axes for open_1d_fits")
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

    XAxis = make_axis(xarr,hdr,wcstype=wcstype,**kwargs)

    return spec,errspec,XAxis,hdr

def make_axis(xarr,hdr,specname=None, wcstype=''):
    """
    Parse parameters from a .fits header into required SpectroscopicAxis
    parameters
    """

    xunits = hdr.get('CUNIT1'+wcstype)
    if hdr.get('ORIGIN') == 'CLASS-Grenoble' and xunits is None:
        # CLASS default
        xunits = 'Hz'

    if hdr.get('REFFREQ'+wcstype):
        reffreq = hdr.get('REFFREQ'+wcstype)
    elif hdr.get('RESTFREQ'+wcstype):
        reffreq = hdr.get('RESTFREQ'+wcstype)
    elif hdr.get('RESTFRQ'+wcstype):
        reffreq = hdr.get('RESTFRQ'+wcstype)
    else:
        reffreq = None

    if hdr.get('CTYPE1'+wcstype):
        xtype = hdr.get('CTYPE1'+wcstype)
    else:
        xtype = 'VLSR'

    XAxis = units.SpectroscopicAxis(xarr,xunits,xtype=xtype,reffreq=reffreq)

    return XAxis

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

    XAxis = make_axis(xarr,hdr,wcstype=wcstype,**kwargs)

    return cube,XAxis,hdr,f
