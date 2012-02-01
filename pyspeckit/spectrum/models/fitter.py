"""
====================
SimpleFitter wrapper 
====================
Adds a variable height (background) component to any model
"""
import numpy
from mpfit import mpfit
from numpy.ma import median
from pyspeckit.spectrum.moments import moments

class SimpleFitter(object):

    def __init__():
        pass

    def moments(self, *args, **kwargs):
        """
        Get the spectral moments from the moments package
        """
        if hasattr(self,'npars') and self.npars == 3:
            return moments(*args,**kwargs)
        elif hasattr(self,'npars') and self.npars > 3:
            m = moments(*args,**kwargs)  
            n_extra_pars = self.npars - len(m) + self.vheight
            return m + ([0] * n_extra_pars)
        else:
            # this may raise an exception itself if npars not defined...
            raise ValueError('Moments not defined for fitter with npars=%i' % self.npars)


def vheightmodel(zeroheightmodel):
    def vhm(xax, *pars,**kwargs):
        """
        Wrapper function vhm to set variable height.
        Parameter order: height, amplitude, shift, width
        """
        vheight=True
        if 'vheight' in kwargs:
            vheight = kwargs.pop('vheight')
        if vheight:
            return zeroheightmodel(xax, *pars[1:],**kwargs) + pars[0]
        else:
            return zeroheightmodel(xax, *pars[1:],**kwargs)
    vhm.__doc__ += zeroheightmodel.__doc__
    return vhm
