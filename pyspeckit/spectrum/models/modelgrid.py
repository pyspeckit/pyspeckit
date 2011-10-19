"""
==========
Model Grid
==========
Fit a line based on parameters output from a grid of models
"""
import numpy as np
from mpfit import mpfit
import matplotlib.cbook as mpcb
import copy
try:
    import scipy.interpolate
    import scipy.ndimage
    scipyOK = True
except ImportError:
    scipyOK=False

def gaussian_line(xax, maxamp, tau, offset, width):
    """
    A Gaussian line function in which the 
    """
    return np.exp(-(xax-offset)**2/(2.0*width**2)) * maxamp * (1.0-np.exp(-1*tau))

def line_params_2D(gridval1, gridval2, griddim1, griddim2, valuegrid):
    """
    Given a 2D grid of modeled line values - the amplitude, e.g. excitation temperature,
    and the optical depth, tau - return the model spectrum

    griddims contains the names of the axes and their values... it should have the same 
    number of entries as gridpars
    """

    if not scipyOK:
        raise ImportError("Scipy could not be imported, therefore interpolation is not available.")

    #matchpt1 = np.argmin( np.abs( gridval1 - griddim1[0,:] ))
    #matchpt2 = np.argmin( np.abs( gridval2 - griddim2[:,0] ))

    return scipy.ndimage.map_coordinates(valuegrid,np.array([[gridval2],[gridval1]]),order=1)

    interpgrid = scipy.interpolate.interp2d( 
            griddim1[ gridval1-5:gridval1+5, gridval2-5:gridval2+5].ravel(),
            griddim2[ gridval1-5:gridval1+5, gridval2-5:gridval2+5].ravel(),
            valuegrid[gridval1-5:gridval1+5, gridval2-5:gridval2+5].ravel())
    
    return interpgrid(gridval1,gridval2)

def line_model_2par(xax, center, width, gridval1, gridval2, griddim1, griddim2, maxampgrid, taugrid,
        linefunction=gaussian_line):
    """
    Returns the spectral line that matches the given x-axis

    xax, center, width must be in the same units!
    """

    maxamp = line_params_2D(gridval1, gridval2, griddim1, griddim2, maxampgrid)
    tau = line_params_2D(gridval1, gridval2, griddim1, griddim2, taugrid)

    return linefunction(xax,maxamp,tau,offset,width)


