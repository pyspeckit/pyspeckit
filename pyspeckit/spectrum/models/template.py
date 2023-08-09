"""
Spectral Template Fitter
========================

A tool to find the optimal shift and scaling for a given template model.

Module API
^^^^^^^^^^
"""
from ..interpolation import interp,_interp
import numpy as np
from . import model

def spectral_template_generator(template_spectrum, xshift_units='km/s', left=0,
                                right=0):
    """
    Given a spectral_template, return a model function with scale and shift as
    free parameters.

    Parameters
    ----------
    template_spectrum: `pyspeckit.spectrum.classes.Spectrum`
        The template spectrum to fit
    xshift_units: str
        The units of the shift parameter
    left/right: float
        The left and right edge parameters used for extrapolating outside the
        template if the template is smaller than the input spectrum.  These
        cannot be NaN.

    Returns
    -------
    spectral_template: function
        The model function that interpolates the template onto the given X-axis
    """
    def spectral_template(xarr, scale, xshift, xshift_units=xshift_units):
        """
        Given a template Spectrum (which should be a Spectrum instance),
        scale & shift it
        """
        shift = xarr.x_to_coord(xshift, xshift_units)
        model = scale * _interp(xarr,
                                template_spectrum.xarr.as_unit(xarr.unit)+shift,
                                template_spectrum.data,
                                left=left,
                                right=right,
                               )

        return model

    return spectral_template

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
    >>> template_fitter = pyspeckit.models.template_fitter(template,
    ...                                                    xshift_units='angstroms')
    >>> dataspec.Registry.add_fitter('template',template_fitter, 2)
    >>> dataspec.specfit(fittype='template',guesses=[1,0])
    >>> print dataspec.specfit.parinfo
    """

    modelfunc = spectral_template_generator(template_spectrum,
                                            xshift_units=xshift_units)

    myclass =  model.SpectralModel(modelfunc, 2, parnames=['scale','shift'],
                                   parlimited=[(True,False),(False,False)],
                                   parlimits=[(0,0), (0,0)],
                                   shortvarnames=('A',r'\Delta x'),
                                   centroid_par='shift',)
    myclass.__name__ = "spectral_template"

    return myclass
