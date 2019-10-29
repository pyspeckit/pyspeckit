from __future__ import print_function
from six.moves import xrange

try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
import numpy as np

from .. import units

def read_galex(fitsfilename, orderselection='obj'):
    """
    Read in a GALEX xg-gsp.fits file 
    Definition of the file is here: 
    `<http://www.galex.caltech.edu/DATA/gr1_docs/GR1_Column_descriptions_for_-xg-gsp.fits.htm>`_

    Parameters
    ----------
    fitsfilename: string 
        The file name (must be a valid FITS file)
    orderselection: {'obj','objs'}
        'obj' selects the NUV 1st and FUV 2nd orders
        'objs' selects the NUV 2nd and FUV 3rd orders
            (higher resolution, lower S/N)

    Example
    -------
    # from http://galex.stsci.edu/GR6/?page=downloadlist&tilenum=5697&type=coaddS
    $ wget http://galex.stsci.edu/data/GR6/pipe/01-vsn/05697-VAR_J000836p545700/g/01-main/0001-img/07-try/VAR_J000836p545700-xg-gsp.fits.gz
    $ gunzip VAR_J000836p545700-xg-gsp.fits.gz

    >>> spblock = read_galex('VAR_J000836p545700-xg-gsp.fits')
    >>> spavg = spblock.average() # average over all spectra of the source (naive)
    >>> spblock.plotter()
    >>> spavg.plotter()
    """
    from .. import classes

    ff = pyfits.open(fitsfilename)
    bintable = ff[1].data

    # wavelength in angstroms
    wavelength = bintable.zero + np.arange(len(bintable.disp)) * bintable.disp
    xaxis = units.SpectroscopicAxis(wavelength,units='angstroms')

    splist = [Spectrum(data=bintable[orderselection][:,ii], error=bintable[orderselection+"err"][:,ii], xarr=xaxis, header=ff[1].header)
            for ii in xrange(bintable[orderselection].shape[1])]

    galexblock = ObsBlock(splist)

    return galexblock
