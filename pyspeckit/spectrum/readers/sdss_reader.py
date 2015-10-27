from __future__ import print_function
from .. import units
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits

def read_sdss(fitsfile, extension=1, data_keyword='flux', error_keyword='ivar', 
        xarr_keyword='loglam', xarr_units='angstrom'):
    """
    Given an SDSS pyfits hdu (read with pyfits.open or astropy.io.fits.open),
    read one of the binary table extensions into a Spectrum

    Parameters
    ----------
    extension : int > 0
        Which BinTableHDU extension to read from (the 0'th is the overall
        header)
    data_keyword : 'data' | 'wdisp' | 'sky' | 'model'
        Which column of the bintable to read in as the spectrum?  Can be any
        keyword actually stored in the BinTableHDU
    error_keyword : 'ivar' | False
        Which keyword to use for the error variable?  Since SDSS stores as
        inverse variance, if 'ivar' is used, will rescale to be 1-sigma errors
        on the data.  If set to False, will be all zeros (so that you can
        read in a model with no errors defined)
    xarr_keyword : 'loglam' 
        Which X-array (dispersion) keyword to use?  If loglam, xarr will be set
        to 10**loglam 
    xarr_units : 'angstroms'
        Units of the X-array.  Unfortunately, this doesn't appear to be stored
        in the header anywhere.
    """
    if isinstance(fitsfile,str):
        pyfits_hdu = pyfits.open(fitsfile)
    else: # no type-checking, just assume it's an HDU
        pyfits_hdu = fitsfile

    hdr = pyfits_hdu[0]._header

    table = pyfits_hdu[extension].data

    data = table[data_keyword]

    if xarr_keyword == 'loglam':
        xarr = 10**table['loglam']
    else:
        xarr = 10**table[xarr_keyword]

    xarr = units.SpectroscopicAxis(xarr,xarr_units,xtype='wavelength')

    if error_keyword == 'ivar':
        error = 1./table[error_keyword]
    elif not error_keyword:
        error = data * 0
    else:
        error = table[error_keyword]

    return data,error,xarr,hdr
