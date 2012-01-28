import pyfits
from .. import units
import numpy.ma as ma
import numpy as np
from . import make_axis
import operator
from pyspeckit.specwarnings import warn

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


def open_1d_pyfits(pyfits_hdu,specnum=0,wcstype='',specaxis="1",errspecnum=None,
        autofix=True, scale_keyword=None, scale_action=operator.div,
        verbose=False, **kwargs):
    """
    This is open_1d_fits but for a pyfits_hdu so you don't necessarily have to
    open a fits file
    """

    # force things that will be treated as strings to be strings
    # this is primarily to avoid problems with variables being passed as unicode
    wcstype = str(wcstype)
    specaxis = str(specaxis)

    hdr = pyfits_hdu._header
    if autofix: 
        for card in hdr.ascardlist():
            try:
                if verbose: card.verify('fix')
                else: card.verify('silentfix')
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
        if hdr.get('WAT0_001') is not None:
            if 'multispec' in hdr.get('WAT0_001'):
                # treat as an Echelle spectrum from  IRAF
                warn("""
This looks like an Echelle spectrum.   You may want to load it
using pyspeckit.wrappers.load_IRAF_multispec.  The file will still
be successfully read if you continue, but the plotting and fitting packages
will run into errors.""")
                return read_echelle(pyfits_hdu)

        if isinstance(specnum,list):
            # allow averaging of multiple spectra (this should be modified
            # - each spectrum should be a Spectrum instance)
            spec = ma.array(data[specnum,:]).mean(axis=0)
        elif isinstance(specnum,int):
            spec = ma.array(data[specnum,:]).squeeze()
        else:
            raise TypeError(
                "Specnum is of wrong type (not a list of integers or an integer)." +
                "  Type: %s" %
                str(type(specnum)))
        if errspecnum is not None:
            # SDSS supplies weights, not errors.    
            if hdr.get('TELESCOP') == 'SDSS 2.5-M':
                errspec = 1. / np.sqrt(ma.array(data[errspecnum,:]).squeeze())
            else:       
                errspec = ma.array(data[errspecnum,:]).squeeze()
        else:
            errspec = spec*0 # set error spectrum to zero if it's not in the data

    elif hdr.get('NAXIS') > 2:
        for ii in xrange(2,hdr.get('NAXIS')):
            # only fail if extra axes have more than one row
            if hdr.get('NAXIS%i' % ii) > 1:
                raise ValueError("Too many axes for open_1d_fits")
        spec = ma.array(data).squeeze()
        if errspecnum is None: 
            errspec = spec*0 # set error spectrum to zero if it's not in the data
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
        if verbose: print "Loading a CLASS .fits spectrum"
        dv = -1*hdr.get('CDELT1')
        if hdr.get('RESTFREQ'):
            v0 = hdr.get('RESTFREQ') + hdr.get('CRVAL1')
        elif hdr.get('RESTF'):
            v0 = hdr.get('RESTF') + hdr.get('CRVAL1')
        else:
            warn("CLASS file does not have RESTF or RESTFREQ")
        p3 = hdr.get('CRPIX1')
    elif hdr.get(str('CD%s_%s%s' % (specaxis,specaxis,wcstype))):
        dv = hdr['CD%s_%s%s' % (specaxis,specaxis,wcstype)]
        v0 = hdr['CRVAL%s%s' % (specaxis,wcstype)]
        p3 = hdr['CRPIX%s%s' % (specaxis,wcstype)]
        hdr.update('CDELT%s' % specaxis,dv)
        if verbose: print "Using the FITS CD matrix.  PIX=%f VAL=%f DELT=%f" % (p3,v0,dv)
    elif hdr.get(str('CDELT%s%s' % (specaxis,wcstype))):
        dv = hdr['CDELT%s%s' % (specaxis,wcstype)]
        v0 = hdr['CRVAL%s%s' % (specaxis,wcstype)]
        p3 = hdr['CRPIX%s%s' % (specaxis,wcstype)]
        if verbose: print "Using the FITS CDELT value.  PIX=%f VAL=%f DELT=%f" % (p3,v0,dv)
    elif len(data.shape) > 1:
        if verbose: print "No CDELT or CD in header.  Assuming 2D input with 1st line representing the spectral axis."
        # try assuming first axis is X axis
        if hdr.get('CUNIT%s%s' % (specaxis,wcstype)):
            xarr = data[0,:]
            spec = data[1,:]
            if data.shape[0] > 2:
                errspec = data[2,:]
        else:
            raise TypeError("Don't know what type of FITS file you've input; "+
                "its header is not FITS compliant and it doesn't look like it "+
                "was written by pyspeckit.")

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

def read_echelle(pyfits_hdu):
    """
    Read an IRAF Echelle spectrum
    
    http://iraf.noao.edu/iraf/ftp/iraf/docs/specwcs.ps.Z
    """

    hdr = pyfits_hdu.header

    WAT1_dict = dict( [s.split('=') for s in hdr.get("WAT1_001").split()] )
    # hdr.get does not preserve whitespace, but whitespace is ESSENTIAL here!
    WAT_string = ""
    ii = 0
    while (hdr.get("WAT2_%03i" % (ii+1))) is not None:
        WAT_string += (hdr.get("WAT2_%03i" % (ii+1))
                + " "*(68-len(hdr.get("WAT2_%03i" % (ii+1)))))
        ii += 1

    WAT_list = WAT_string.split("spec")

    if "multi" in WAT_list[0]:
        WAT_list.pop(0)
    if ' ' in WAT_list:
        WAT_list.remove(' ')
    if '' in WAT_list:
        WAT_list.remove('')

    x_axes = []

    specaxdict = dict( [ (int(s.split("=")[0]),s.split("=")[1]) for s in WAT_list ] )
    for specnum, axstring in specaxdict.iteritems():
        axsplit = axstring.replace('"','').split()
        if specnum != int(axsplit[0]):
            raise ValueError("Mismatch in IRAF Echelle specification")
        num,beam,dtype,crval,cdelt,naxis,z,aplow,aphigh = axsplit[:9]
        if len(axsplit) > 9:
            functions = axsplit[9:]

        if int(dtype) == 0:
            xax = ( float(crval) + float(cdelt) * (np.arange(int(naxis)) + 1) ) / (1.+float(z))

        headerkws = {'CRPIX1':1, 'CRVAL1':crval, 'CDELT1':cdelt,
                'NAXIS1':naxis, 'NAXIS':1, 'REDSHIFT':z,
                'CTYPE1':'wavelength', 'CUNIT1':WAT1_dict['units']}
        cards = [pyfits.Card(k,v) for (k,v) in headerkws.iteritems()]
        header = pyfits.Header(cards)

        xarr = make_axis(xax,header)
        x_axes.append(xarr)
    
    return pyfits_hdu.data, pyfits_hdu.data*0, units.EchelleAxes(x_axes), hdr
