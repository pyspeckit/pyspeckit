import pyfits
import pyspeckit
from pyspeckit.spectrum.readers.fits_reader import read_echelle


def load_IRAF_multispec(fitsfilename):
    """
    Load an IRAF multispec file as a list of spectra
    (a list of spectra can be passed to `pyspeckit.ObsBlock` if they are all of
    the same wavelength regime, but different objects or different spectra of
    the same object.  They can be passed to `pyspeckit.Spectra` if, e.g., you
    have an Echelle spectrum 
    """

    fitsfile = pyfits.open(fitsfilename)

    data, err, xax, header = read_echelle(fitsfile[0])

    speclist = [pyspeckit.Spectrum(data=data[ii,:], 
            error=err[ii,:],
            xarr=xax[ii,:], 
            header=header) 
        for ii in xrange(data.shape[0])]

    return speclist
