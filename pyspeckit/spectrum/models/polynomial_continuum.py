"""
====================
Polynomial Continuum
====================

A simple polynomial (or linear, lines are polynomials) continuum model
You can access this same model through baseline, but it's here in case the
continuum is what you're really trying to measure.

"""
import model
import numpy 

def polymodel(x, pars,  return_components=False):
    """
    Wrapper of :meth:`numpy.polyval`.  Extra parameter - return_components - is
    ignored, but required.
    """
    __doc__ += numpy.polyval.__doc__
    return numpy.polyval(x, pars)


def poly_fitter(order=1, multisingle='multi'):
    """
    Generator for polynomial fitter class
    """

    myclass =  model.SpectralModel(polymodel, order,
            parnames=['coeff%i' % ii for ii in xrange(order+1)], 
            parlimited=[(False,False) for ii in xrange(order+1)], 
            parlimits=[(0,0) for ii in xrange(order+1)], 
            shortvarnames=('C' for ii in xrange(order+1)),
            multisingle=multisingle,
            )
    myclass.__name__ = "polymodel"
    
    return myclass
