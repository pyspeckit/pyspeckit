"""
====================
Polynomial Continuum
====================

A simple polynomial (or linear, lines are polynomials) continuum model
You can access this same model through baseline, but it's here in case the
continuum is what you're really trying to measure.

Module API
^^^^^^^^^^
"""
from six.moves import xrange
from . import model
import numpy 

def polymodel(x, *pars, **kwargs):
    """
    Wrapper of :meth:`numpy.polyval`.  Extra parameter - return_components - is
    ignored, but required.
    *args = polynomial parameters
    **kwargs = just to catch extra junk; not passed
    """
    # polyval and astropy quantity are incompatible
    if hasattr(x, 'value'):
        x = x.value

    return numpy.polyval(pars, x)

polymodel.__doc__ += numpy.polyval.__doc__

def poly_fitter(order=1):
    """
    Generator for polynomial fitter class
    """

    myclass =  model.SpectralModel(polymodel, order+1,
            parnames=['coeff%i' % ii for ii in xrange(order+1)], 
            parlimited=[(False,False) for ii in xrange(order+1)], 
            parlimits=[(0,0) for ii in xrange(order+1)], 
            shortvarnames=['C%i' % ii for ii in xrange(order+1)]
            )
    myclass.__name__ = "polymodel"
    
    return myclass
