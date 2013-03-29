from ..interpolation import interp,_interp
import numpy as np
import model

def spectral_template_generator(template_spectrum, xshift_units='km/s'):
    """
    Given a spectral_template, return a model function with scale and shift as
    free parameters.
    """
    def spectral_template(xarr,  scale, xshift, xshift_units=xshift_units):
        """
        Given a template Spectrum (which should be a Spectrum instance),
        scale & shift it
        """
        shift = xarr.x_to_units(xshift, xshift_units)
        model = _interp(xarr, template_spectrum.xarr.as_unit(xarr.units)+shift,
                template_spectrum.data, left=np.nan, right=np.nan) * scale
        
        return model

def template_fitter(template_spectrum, xshift_units='km/s'):
    """
    Generator for Spectral Template fitter class

    Parameters
    ----------
    template_spectrum : pyspeckit.Spectrum
        A valid spectrum to be scaled and shifted to match the input
    xshift_units : str in pyspeckit.units.unit_type_dict
        The units of the shift to fit.  If you're using a velocity unit, make
        sure there's a reference X-unit for both the template spectrum and the
        input spectrum.

    Examples
    --------
    >>> template = pyspeckit.Spectrum("template_spectrum.fits")
    >>> dataspec = pyspeckit.Spectrum("DataSpectrum.fits")
    >>> template_fitter = pyspeckit.models.template_fitter(template,xshift_units='angstroms')
    >>> dataspec.Registry.add_fitter('template',template_fitter,2,multisingle='multi')
    >>> dataspec.specfit(fittype='template',guesses=[1,0])
    >>> print dataspec.specfit.parinfo
    """

    modelfunc = spectral_template_generator(template_spectrum, xshift_units=xshift_units)

    myclass =  model.SpectralModel(modelfunc, 2,
            parnames=['scale','shift'], 
            parlimited=[(True,False),(False,False)], 
            parlimits=[(0,0), (0,0)],
            shortvarnames=('A',r'\Delta x'),
            centroid_par='shift',
            )
    myclass.__name__ = "spectral_template"
    
    return myclass
