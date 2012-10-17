"""
====================
Voigt Profile Fitter
====================
"""
import model
import numpy as np
try:
    import scipy.special
    scipyOK = True
except ImportError:
    scipyOK = False

def voigt(xarr,amp,xcen,sigma,gamma,normalized=False):
    """
    Normalized Voigt profile

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
    """

    if scipyOK:
        z = ((xarr-xcen) + 1j*gamma) / (sigma * np.sqrt(2))
        V = amp * np.real(scipy.special.wofz(z)) 
        if normalized:
            return V / (sigma*np.sqrt(2*np.pi))
        else:
            return V
        #tmp = 1.0/scipy.special.wofz(numpy.zeros((len(xarr))) \
        #      +1j*numpy.sqrt(numpy.log(2.0))*Lfwhm).real
        #tmp = tmp*amp* \
        #      scipy.special.wofz(2*numpy.sqrt(numpy.log(2.0))*(xarr-xcen)/Gfwhm+1j* \
        #      numpy.sqrt(numpy.log(2.0))*Lfwhm).real
        #return tmp
    else:
        raise ImportError("Couldn't import scipy, therefore cannot do voigt profile stuff")

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

def voigt_fitter(multisingle='multi'):
    """
    Generator for voigt fitter class
    """

    myclass =  model.SpectralModel(voigt, 4,
            parnames=['amplitude','shift','gwidth','lwidth'], 
            parlimited=[(False,False),(False,False),(True,False),(True,False)], 
            parlimits=[(0,0), (0,0), (0,0), (0,0)],
            shortvarnames=('A',r'\Delta x',r'\sigma_G',r'\sigma_L'),
            multisingle=multisingle,
            centroid_par='shift',
            fwhm_func=voigt_fwhm,
            fwhm_pars=['gwidth','lwidth'],
            )
    myclass.__name__ = "voigt"
    
    return myclass
