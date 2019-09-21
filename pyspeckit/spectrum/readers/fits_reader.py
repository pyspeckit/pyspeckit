from __future__ import print_function
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
import numpy.ma as ma
import numpy as np
from six import operator

from .. import units
from . import make_axis
from ...specwarnings import warn

def open_1d_fits(filename, hdu=0, **kwargs):
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

    # try to open as an HDU...
    if all((hasattr(filename,k) for k in ('data','header','_header'))):
        f = [filename]
        hdu = 0
    elif isinstance(filename,pyfits.HDUList):
        f = filename
    else:
        f = pyfits.open(filename, ignore_missing_end=True)

    return open_1d_pyfits(f[hdu],**kwargs)


def open_1d_pyfits(pyfits_hdu, specnum=0, wcstype='', specaxis="1",
                   errspecnum=None, autofix=True, scale_keyword=None,
                   scale_action=operator.truediv, verbose=False, apnum=0,
                   **kwargs):
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
        for card in hdr.cards:
            try:
                if verbose:
                    card.verify('fix')
                else:
                    card.verify('silentfix')
            except pyfits.VerifyError:
                del hdr[card.keyword]

    data = pyfits_hdu.data

    with np.errstate(invalid='ignore'):
        # silently turn signalling nans into quiet nans
        data[np.isnan(data)] = np.nan

    # search for the correct axis (may be 1 or 3, unlikely to be 2 or others)
    # 1 = 1D spectrum
    # 3 = "3D" spectrum with a single x,y point (e.g., JCMT smurf/makecube)
    if hdr.get('NAXIS') > 1:
        for ii in range(1,hdr.get('NAXIS')+1):
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
        if hdr.get('BANDID2'):
            # this is an IRAF .ms.fits file with a 'background' in the 3rd dimension
            spec = ma.array(data[specnum,apnum,:]).squeeze()
        else:
            for ii in range(1,hdr.get('NAXIS')+1):
                # only fail if extra axes have more than one row
                if hdr.get('NAXIS%i' % ii) > 1 and (ii != int(specaxis)):
                    raise ValueError("Too many axes for open_1d_fits")
            spec = ma.array(data).squeeze()
        if errspecnum is None:
            errspec = spec*0 # set error spectrum to zero if it's not in the data
    else:
        spec = ma.array(data).squeeze()
        if errspecnum is None: errspec = spec*0 # set error spectrum to zero if it's not in the data

    if scale_keyword is not None:
        try:
            print(("Found SCALE keyword %s.  Using %s to scale it" % (scale_keyword,scale_action)))
            scaleval = hdr[scale_keyword]
            spec = scale_action(spec,scaleval)
            errspec = scale_action(errspec,scaleval)
        except (ValueError, KeyError) as e:
            pass

    xarr = None
    if hdr.get('ORIGIN') == 'CLASS-Grenoble':
        # Use the CLASS FITS definition (which is non-standard)
        # http://iram.fr/IRAMFR/GILDAS/doc/html/class-html/node84.html
        # F(n) = RESTFREQ + CRVALi + ( n - CRPIXi ) * CDELTi
        if verbose: print("Loading a CLASS .fits spectrum")
        # there is no reason to assume CDELT is wrong dv = -1*hdr.get('CDELT1')
        dv = hdr.get('CDELT1')
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
        hdr['CDELT%s' % specaxis] = dv
        if verbose: print(("Using the FITS CD matrix.  PIX=%f VAL=%f DELT=%f" % (p3,v0,dv)))
    elif hdr.get(str('CDELT%s%s' % (specaxis,wcstype))):
        if hdr.get(str('PC{0}_{0}'.format(specaxis))):
            dv = hdr['CDELT%s%s' % (specaxis,wcstype)] * hdr.get(str('PC{0}_{0}'.format(specaxis)))
        else:
            dv = hdr['CDELT%s%s' % (specaxis,wcstype)]
        v0 = hdr['CRVAL%s%s' % (specaxis,wcstype)]
        p3 = hdr['CRPIX%s%s' % (specaxis,wcstype)]
        if verbose: print(("Using the FITS CDELT value.  PIX=%f VAL=%f DELT=%f" % (p3,v0,dv)))
    elif len(data.shape) > 1:
        if verbose: print("No CDELT or CD in header.  Assuming 2D input with 1st line representing the spectral axis.")
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

    # need to do something with this...
    restfreq = hdr.get('RESTFREQ')
    if restfreq is None: restfreq= hdr.get('RESTFRQ')

    XAxis = make_axis(xarr,hdr,wcstype=wcstype,specaxis=specaxis,**kwargs)

    return spec,errspec,XAxis,hdr


def _get_WATS(hdr):
    """
    Get the WATS from an IRAF Echelle header
    """
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

    specaxdict = dict( [ (int(s.split("=")[0]),s.split("=")[1]) for s in WAT_list ] )

    return WAT1_dict,specaxdict

def make_linear_axis(hdr, axsplit, WAT1_dict):
    num,beam,dtype,crval,cdelt,naxis,z,aplow,aphigh = axsplit[:9]
    # this is a hack for cropped spectra...
    #print "header naxis: %i, WAT naxis: %i" % (hdr['NAXIS1'], int(naxis))
    if hdr['NAXIS1'] != int(naxis):
        naxis = hdr['NAXIS1']
        crpix = hdr.get('CRPIX1')
        warn("Treating as cropped echelle spectrum.")
    else:
        crpix = 0

    if len(axsplit) > 9:
        functions = axsplit[9:]
        warn("Found but did not use functions %s" % str(functions))

    if int(dtype) == 0:
        # Linear dispersion (eq 2, p.5 from Valdez, linked above)
        xax = ((float(crval) + float(cdelt) * (np.arange(int(naxis)) + 1 -
                                               crpix)) / (1.+float(z)))
    else: raise ValueError("Unrecognized LINEAR dispersion in IRAF Echelle specification")

    headerkws = {'CRPIX1':1, 'CRVAL1':crval, 'CDELT1':cdelt,'NAXIS1':naxis,
             'NAXIS':1, 'REDSHIFT':z, 'CTYPE1':'wavelength',
             'CUNIT1':WAT1_dict['units']}

    return xax, naxis, headerkws

def make_multispec_axis(hdr, axsplit, WAT1_dict):
    num,beam,dtype,crval,cdelt,naxis,z,aplow,aphigh = axsplit[:9]
    # this is a hack for cropped spectra...
    #print "header naxis: %i, WAT naxis: %i" % (hdr['NAXIS1'], int(naxis))
    if hdr['NAXIS1'] != int(naxis):
        crpix = int(naxis) - hdr['NAXIS']
        naxis = hdr['NAXIS1']
        warn("Treating as cropped echelle spectrum.")
    else:
        crpix = 0

    if len(axsplit) > 9:
        functions = axsplit[9:]
        warn("Found but did not use functions %s" % str(functions))

    if int(dtype) == 0:
        # Linear dispersion (eq 11, p.9 from Valdez, linked above)
        xax = (float(crval) + float(cdelt) * (np.arange(int(naxis)))) / (1.+float(z))
    elif int(dtype) == 1:
        # Log-linear dispersion (eq 12, p.9 from Valdez, linked above)
        xax = 10.**(float(crval) + float(cdelt) * (np.arange(int(naxis)))) / (1.+float(z))
    # elif int(dtype) == 2:
        # Non-linear dispersion
    # elif int(dtype) == -1:
        # Data is not dispersion coords
    else: raise ValueError("Unrecognized MULTISPE dispersion in IRAF Echelle specification")

    headerkws = {'CRPIX1':1, 'CRVAL1':crval, 'CDELT1':cdelt, 'NAXIS1':naxis,
                 'NAXIS':1, 'REDSHIFT':z, 'CTYPE1':'wavelength',
                 'CUNIT1':WAT1_dict['units']}

    return xax, naxis, headerkws

def read_echelle(pyfits_hdu):
    """
    Read an IRAF Echelle spectrum

    http://iraf.noao.edu/iraf/ftp/iraf/docs/specwcs.ps.Z
    """

    hdr = pyfits_hdu.header

    WAT1_dict, specaxdict = _get_WATS(hdr)

    x_axes = []

    for specnum, axstring in specaxdict.items():
        axsplit = axstring.replace('"','').split()
        if specnum != int(axsplit[0]):
            raise ValueError("Mismatch in IRAF Echelle specification")
        num,beam,dtype,crval,cdelt,naxis,z,aplow,aphigh = axsplit[:9]

        if hdr['CTYPE1'] == 'LINEAR':
            xax,naxis,headerkws = make_linear_axis(hdr, axsplit, WAT1_dict)
        elif hdr['CTYPE1'] == 'MULTISPE':
            xax,naxis,headerkws = make_multispec_axis(hdr, axsplit, WAT1_dict)

        cards = [pyfits.Card(k,v) for (k,v) in headerkws.items()]
        header = pyfits.Header(cards)

        xarr = make_axis(xax,header)
        x_axes.append(xarr)

    data = pyfits_hdu.data
    with np.errstate(invalid='ignore'):
        data[np.isnan(data)] = np.nan

    return data, np.zeros_like(data), units.EchelleAxes(x_axes), hdr
