"""
====================
Voigt Profile Fitter
====================
"""
import model
import numpy
try:
    import scipy.special
    scipyOK = True
except ImportError:
    scipyOK = False

def voigt(xarr,amp,xcen,Gfwhm,Lfwhm):
    """
    voigt profile

    V(x,sig,gam) = Re(w(z))/(sig*sqrt(2*pi))
    z = (x+i*gam)/(sig*sqrt(2))

    Converted from 
    http://mail.scipy.org/pipermail/scipy-user/2011-January/028327.html
    """

    if scipyOK:
        tmp = 1.0/scipy.special.wofz(numpy.zeros((len(xarr))) \
              +1j*numpy.sqrt(numpy.log(2.0))*Lfwhm).real
        tmp = tmp*amp* \
              scipy.special.wofz(2*numpy.sqrt(numpy.log(2.0))*(xarr-xcen)/Gfwhm+1j* \
              numpy.sqrt(numpy.log(2.0))*Lfwhm).real
        return tmp
    else:
        raise ImportError("Couldn't import scipy, therefore cannot do voigt profile stuff")

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
            )
    myclass.__name__ = "voigt"
    
    return myclass
