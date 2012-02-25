"""
ALFAFA "source" .sav file
"""
import idlsave
import pyfits
import pyspeckit
import numpy as np

def read_alfalfa(filename, sourcenumber=0):
    """
    Create an Observation Block class for a single source in an ALFALFA
    'source' IDL save file
    """
    
    savfile = idlsave.read(filename)
    src = savfile.src[sourcenumber]

    header = pyfits.Header()

    splist = []
    for spectra in src.spectra:
        for par in spectra.dtype.names:
            try:
                len(spectra[par])
            except TypeError:
                header.update(par[:8],spectra[par])

        # assume ALFALFA spectra in Kelvin
        header.update('BUNIT','K')

        xarr = pyspeckit.spectrum.units.SpectroscopicAxis(spectra.velarr,
                refX=header['RESTFRQ'], refX_units='MHz', unit='km/s')

        data = np.ma.masked_where(np.isnan(spectra.spec), spectra.spec)

        sp = pyspeckit.Spectrum(xarr=xarr, data=data, header=header)

        # the Source has a baseline presubtracted (I think)
        sp.baseline.baselinepars = spectra.baseline[::-1]
        sp.baseline.subtracted = True
        sp.baseline.order = len(spectra.baseline)
        sp.baseline.basespec = np.poly1d(sp.baseline.baselinepars)(np.arange(xarr.shape[0]))

        # There are multiple components in each Spectrum, but I think they are not indepedent
        sp.specfit.fittype = 'gaussian'
        sp.specfit.fitter = sp.specfit.Registry.multifitters['gaussian']
        modelpars = zip(spectra['STON'],spectra['VCEN'],spectra['WIDTH']) 
        modelerrs = zip(spectra['STON'],spectra['VCENERR_STAT'],spectra['WIDTHERR']) 
        sp.specfit.modelpars = modelpars[0] # only use the first fit
        sp.specfit.fitter.mpp = modelpars[0]
        sp.specfit.modelerrs = modelerrs[0]
        sp.specfit.fitter.mpperr = modelerrs[0]
        sp.specfit._full_model()

        splist.append(sp)

    return pyspeckit.ObsBlock(splist)
