"""
====================
Power-law Continuum
====================

A simple power-law continuum model
You can access this same model through baseline, but it's here in case the
continuum is what you're really trying to measure.

"""
from . import model
import numpy 

def powerlaw(x, scale, power, return_components=False):
    """
    Defines a power law

    Returns
    -------
    scale * x**power
    """
    return scale*x**power


def powerlaw_fitter(order=1):
    """
    Generator for powerlaw fitter class
    """

    myclass =  model.SpectralModel(powerlaw, 2,
            parnames=['scale','power'], 
            parlimited=[(False,False),(False,False)], 
            parlimits=[(0,0),(0,0)], 
            shortvarnames=('S','P')
            )
    myclass.__name__ = "powerlaw"
    
    return myclass
