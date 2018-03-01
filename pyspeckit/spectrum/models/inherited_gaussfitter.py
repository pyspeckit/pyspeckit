"""
===============
Gaussian Fitter
===============

The simplest and most useful model.

Until 12/23/2011, gaussian fitting used the complicated and somewhat bloated
gaussfitter.py code.  Now, this is a great example of how to make your own
model!  Just make a function like gaussian and plug it into the SpectralModel
class.

Module API
^^^^^^^^^^
"""
import numpy as np
import warnings
from . import model
from . import fitter
from ...specwarnings import PyspeckitWarning

def gaussian(xarr, amplitude, dx, width, return_components=False, normalized=False,
             return_hyperfine_components=False):
    """
    Returns a 1-dimensional gaussian of form
    A*np.exp(-(x-dx)**2/(2*w**2))

    Area is sqrt(2*pi*sigma^2)*amplitude - i.e., this is NOT a normalized
    gaussian, unless normalized=True in which case A = Area
    
    Parameters
    ----------
    xarr : np.ndarray
        array of x values
    amplitude : float
        Amplitude of the Gaussian, i.e. its peak value, unless
        normalized=True then A is the area of the gaussian
    dx : float
        Center or "shift" of the gaussian
    width : float
        Width of the gaussian (sigma)
    return_components : bool
        dummy variable; return_components does nothing but is required by all
        fitters
    return_hyperfine_components : bool
        dummy variable; does nothing but is required by all
        fitters
    normalized : bool
        Return a normalized Gaussian?
    """
    if width == 0:
        return np.nan
    elif width < 0:
        warnings.warn("Negative width in Gaussian: {0}.".format(width),
                      PyspeckitWarning)

    xarr = np.array(xarr) # make sure xarr is no longer a spectroscopic axis
    model = amplitude*np.exp(-(xarr-dx)**2/(2.0*width**2))
    if normalized:
        return model / (np.sqrt(2*np.pi) * width**2)
    else:
        return model

def gaussian_fwhm(sigma):
    return np.sqrt(8*np.log(2)) * sigma

def gaussian_integral(amplitude, sigma):
    """ Integral of a Gaussian """
    return amplitude * np.sqrt(2*np.pi*sigma**2)

def _integral_modelpars(modelpars=None):
    """ light wrapper to match requirements for model.analytic_integral """
    amplitude = modelpars[0]
    sigma = modelpars[2]
    return gaussian_integral(amplitude,sigma)

def gaussian_fitter():
    """
    Generator for Gaussian fitter class
    """

    myclass = model.SpectralModel(gaussian, 3,
            parnames=['amplitude','shift','width'],
            parlimited=[(False,False),(False,False),(True,False)],
            parlimits=[(0,0), (0,0), (0,0)],
            shortvarnames=('A',r'\Delta x',r'\sigma'),
            centroid_par='shift',
            fwhm_func=gaussian_fwhm,
            fwhm_pars=['width'],
            integral_func=_integral_modelpars,
            )
    myclass.__name__ = "gaussian"
    
    return myclass

def gaussian_vheight_fitter():
    """
    Generator for Gaussian fitter class
    """

    vhg = fitter.vheightmodel(gaussian)
    myclass = model.SpectralModel(vhg, 4,
            parnames=['height','amplitude','shift','width'],
            parlimited=[(False,False),(False,False),(False,False),(True,False)],
            parlimits=[(0,0),(0,0), (0,0), (0,0)],
            shortvarnames=('B','A',r'\Delta x',r'\sigma'),
            centroid_par='shift',
            fwhm_func=gaussian_fwhm,
            fwhm_pars=['width'],
            )
    myclass.__name__ = "vheightgaussian"
    
    return myclass
