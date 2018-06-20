"""
====================
Voigt Profile Fitter
====================

Module API
^^^^^^^^^^
"""
from . import model
import numpy as np
from ...spectrum.moments import moments
import types
try:
    import scipy.special
    scipyOK = True
except ImportError:
    scipyOK = False


def voigt(xarr, amp, xcen, sigma, gamma, normalized=False):
    """
    Normalized Voigt profile

    z = (x+i*gam)/(sig*sqrt(2))
    V(x,sig,gam) = Re(w(z))/(sig*sqrt(2*pi))

    The area of V in this definition is 1.
    If normalized=False, then you can divide the integral of V by
    sigma*sqrt(2*pi) to get the area.

    Original implementation converted from
    http://mail.scipy.org/pipermail/scipy-user/2011-January/028327.html
    (had an incorrect normalization and strange treatment of the input
    parameters)

    Modified implementation taken from wikipedia, using the definition.
    http://en.wikipedia.org/wiki/Voigt_profile

    Parameters
    ----------
    xarr : np.ndarray
        The X values over which to compute the Voigt profile
    amp : float
        Amplitude of the voigt profile
        if normalized = True, amp is the AREA
    xcen : float
        The X-offset of the profile
    sigma : float
        The width / sigma parameter of the Gaussian distribution
    gamma : float
        The width / shape parameter of the Lorentzian distribution
    normalized : bool
        Determines whether "amp" refers to the area or the peak
        of the voigt profile
    """

    if scipyOK:
        z = ((xarr.value-xcen) + 1j*gamma) / (sigma * np.sqrt(2))
        V = amp * np.real(scipy.special.wofz(z))
        if normalized:
            return V / (sigma*np.sqrt(2*np.pi))
        else:
            return V
    else:
        raise ImportError("Couldn't import scipy, therefore cannot do "
                          "voigt profile stuff")


def voigt_fwhm(sigma, gamma):
    """
    Approximation to the Voigt FWHM from wikipedia

    http://en.wikipedia.org/wiki/Voigt_profile

    Parameters
    ----------
    sigma : float
        The width / sigma parameter of the Gaussian distribution
    gamma : float
        The width / shape parameter of the Lorentzian distribution
    """
    return 0.5346 * 2 * gamma + np.sqrt(0.2166*(2*gamma)**2 + sigma**2*8*np.log(2))


def voigt_moments(self, *args, **kwargs):
    """
    Get the spectral moments from the moments package.  Use the gaussian width
    for the lorentzian width (not a great guess!)
    """
    m = moments(*args,**kwargs)
    return list(m) + [m[-1]]


def voigt_fitter():
    """
    Generator for voigt fitter class
    """

    myclass = model.SpectralModel(voigt, 4,
                                  parnames=['amplitude', 'shift', 'gwidth',
                                            'lwidth'],
                                  parlimited=[(False, False), (False, False),
                                              (True, False), (True, False)],
                                  parlimits=[(0, 0),  (0, 0),  (0, 0),  (0, 0)],
                                  shortvarnames=('A', r'\Delta x',
                                                 r'\sigma_G', r'\sigma_L'),
                                  centroid_par='shift',
                                  fwhm_func=voigt_fwhm,
                                  fwhm_pars=['gwidth','lwidth'],
                                  guess_types=['amplitude', 'center', 'width',
                                               'width'],
            )
    myclass.__name__ = "voigt"
    try:
        myclass.moments = types.MethodType(voigt_moments, myclass,
                                           myclass.__class__)
    except TypeError: # indicates py3 is being used
        # http://stackoverflow.com/questions/10729909/convert-builtin-function-type-to-method-type-in-python-3?lq=1
        myclass.moments = voigt_moments.__get__(myclass)

    return myclass
