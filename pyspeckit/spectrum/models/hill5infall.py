"""
===========================
Hill5 analytic infall model 
===========================
.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>

Code translated from:
https://bitbucket.org/devries/analytic_infall/overview

Original source:
http://adsabs.harvard.edu/abs/2005ApJ...620..800D
"""
import numpy as np
import model

def hill5_model( xarr, tau, v_lsr,  v_infall,  sigma,  tpeak, TBG=2.73):
    """
    The rest of this needs to be translated from C
    """
  
    vf = v_lsr+v_infall
    vr = v_lsr-v_infall
  
    velocity_array = np.array(xarr.as_unit('km/s'))
    tauf = tau*np.exp(-((velocity_array-vf)/sigma)**2 / 2.0 )
    taur = tau*np.exp(-((velocity_array-vr)/sigma)**2 / 2.0 )

    subf = ( (1-np.exp(-tauf))/tauf * (tauf>1e-4) + (tauf < 1e-4) )
    subr = ( (1-np.exp(-taur))/taur * (taur>1e-4) + (taur < 1e-4) )
    subf[subf != subf] = 1.0
    subr[subr != subr] = 1.0
  
    frequency = np.array(xarr.as_unit('Hz'))
    hill_array =(jfunc(tpeak,frequency)-jfunc(TBG,frequency))*(subf-np.exp(-tauf)*subr)

    if np.isnan(hill_array).any():
        raise ValueError("Hill5 model has a NAN")

    return hill_array

def jfunc(t, nu):
    """
    t- kelvin
    nu - Hz?
    """


    H=6.6260755e-27
    K=1.380658e-16
    NSIG=2.0

    # I think this is only necessary for C
    #if(nu<1.0e-6) return t
    
    to = H*nu/K
    return to/(np.exp(to/t)-1.0)

hill5_fitter = model.SpectralModel(hill5_model, 5,
        parnames=['tau', 'v_lsr',  'v_infall',  'sigma', 'tpeak'], 
        parlimited=[(True,False),(False,False),(True,False),(True,False), (True,False)], 
        parlimits=[(0,0), (0,0), (0,0), (0,0), (0,0)],
        shortvarnames=("\\tau","v_{lsr}","v_{infall}","\\sigma","T_{peak}"), # specify the parameter names (TeX is OK)
        fitunits='Hz' )
