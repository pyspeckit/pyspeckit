"""
===========================
Formaldehyde mm-line fitter
===========================

This is a formaldehyde 3_03-2_02 / 3_22-221 and 3_03-2_02/3_21-2_20 fitter.
It is based entirely on RADEX models.

This is the EWR fork of the fitter in pyspeckit.

Module API
^^^^^^^^^^
"""
from __future__ import print_function
import numpy as np
from . import hyperfine
from . import fitter,model#,modelgrid
try: # for model grid reading
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
try:
    import scipy.interpolate
    import scipy.ndimage
    scipyOK = True
except ImportError:
    scipyOK=False
from six.moves import xrange

# h2co_mm_vtau = hyperfine.hyperfinemodel(line_names, voff_lines_dict,
#         freq_dict, line_strength_dict, relative_strength_total_degeneracy)
# h2co_mm_vtau_fitter = h2co_mm_vtau.fitter
# h2co_mm_vtau_vheight_fitter = h2co_mm_vtau.vheight_fitter



line_names = ['threeohthree','threetwotwo','threetwoone']

# http://adsabs.harvard.edu/abs/1971ApJ...169..429T has the most accurate freqs
# http://adsabs.harvard.edu/abs/1972ApJ...174..463T [twotwo]
central_freq_dict = {
        'threeohthree': 218.222192e9,
        'threetwotwo': 218.475632e9,
        'threetwoone': 218.760066e9,
    }
line_strength_dict={
        'threeohthree': 1.,
        'threetwotwo': 1.,
        'threetwoone': 1.,
        }
relative_strength_total_degeneracy={
        'threeohthree': 1.,
        'threetwotwo': 1.,
        'threetwoone': 1.,
        }
freq_dict = central_freq_dict
aval_dict = {
        'threeohthree': 2.818e-4,
        'threetwotwo': 1.571e-4,
        'threetwoone': 1.577e-4,
    }

voff_lines_dict = {
        'threeohthree': 0.,
        'threetwotwo': 0.,
        'threetwoone': 0.,
        }


formaldehyde_mm_vtau = hyperfine.hyperfinemodel(line_names, voff_lines_dict,
        freq_dict, line_strength_dict, relative_strength_total_degeneracy)
formaldehyde_mm_vtau_fitter = formaldehyde_mm_vtau.fitter
formaldehyde_mm_vtau_vheight_fitter = formaldehyde_mm_vtau.vheight_fitter

def h2co_mm_radex(xarr,
        Temperature=25,
        logColumn=13,
        logDensity=4,
        xoff_v=0.0,
        width=1.0,
        grid_vwidth=1.0,
        gridbundle = None,
        debug=False,
        verbose=False,
        **kwargs):

    """
    Use a grid of RADEX-computed models to make a model line spectrum

    The RADEX models have to be available somewhere.
    OR they can be passed as arrays.  If as arrays, the form should be:
    texgrid = ((minfreq1,maxfreq1,texgrid1),(minfreq2,maxfreq2,texgrid2))

    xarr must be a SpectroscopicAxis instance
    xoff_v, width are both in km/s


    Parameters
    ----------
    grid_vwidth : float
        the velocity assumed when computing the grid in km/s
        this is important because tau = modeltau / width (see, e.g.,
        Draine 2011 textbook pgs 219-230)
    density : float
        Density!
    """

    # Convert X-units to frequency in GHz
    xarr = xarr.as_unit('Hz', quiet=True)
    Tex303,Tex322,Tex321,tau303,tau322,tau321 = gridbundle


    # if this gets too far different from 1, we are gonna have a Bad Time.
    scalefac = grid_vwidth/width

    tex = (Tex303(logColumn,logDensity,Temperature),
           Tex322(logColumn,logDensity,Temperature),
           Tex321(logColumn,logDensity,Temperature))
    tau =  (tau303(logColumn,logDensity,Temperature)*scalefac,
            tau322(logColumn,logDensity,Temperature)*scalefac,
            tau321(logColumn,logDensity,Temperature)*scalefac)
    if np.any(np.isnan(tex)) or np.any(np.isnan(tau)):
        raise ValueError("Invalid column/density")

    if verbose:
        for ta,tk in zip(tau,tex):
            print("density %20.12g temperature %20.12g column %20.12g: tau %20.12g tex %20.12g" % (logDensity, Temperature, logColumn, ta, tk))

    if debug:
        import pdb; pdb.set_trace()

# here there be physics
    ckms = 2.99792458e5
    freq_dict = {
        '303': 218.222192e9,
        '322': 218.475632e9,
        '321': 218.760066e9,
        }
    Tbg = 2.73 #because it totally is


    nu0 = np.array([ 218.222192e9, 218.475632e9,218.760066e9])
    nuwidth = [width/ckms*nu for nu in nu0]
    nuoff = [xoff_v/ckms*nu for nu in nu0]
    minfreq = nu0/1e9 - 0.25
    maxfreq = nu0/1e9 + 0.25
#    spec2 = np.zeros(len(xarr))
#    for ii in range(len(nu0)):
#        taunu = tau[ii]*np.exp(-(xarr+nuoff[ii]-nu0[ii])**2/(2.0*nuwidth[ii]**2))
#        spec2 = spec2 + (1-np.exp(-taunu))*tex[ii] + Tbg*(np.exp(-taunu)-1)  #second term assumes an ON-OFF

    spec = np.sum([
            (formaldehyde_mm_vtau(xarr, Tex=float(tex[ii]), tau=float(tau[ii]),
                                  xoff_v=xoff_v, width=width, **kwargs)
             * (xarr.as_unit('GHz')>minfreq[ii]) * (xarr.as_unit('GHz')<maxfreq[ii])) for ii in xrange(len(tex))],
                  axis=0)
#    import pdb
#    pdb.set_trace()


    return spec


def formaldehyde_mm(xarr, amp=1.0, xoff_v=0.0, width=1.0,
        return_components=False ):
    """
    Generate a model Formaldehyde spectrum based on simple gaussian parameters

    the "amplitude" is an essentially arbitrary parameter; we therefore define it to be Tex given tau=0.01 when
    passing to the fitter
    The final spectrum is then rescaled to that value

    The components are independent, but with offsets set by frequency... in principle.

    """

    mdl = formaldehyde_vtau(xarr, Tex=amp*0.01, tau=0.01, xoff_v=xoff_v,
            width=width,
            return_components=return_components)
    if return_components:
        mdlpeak = np.abs(mdl).squeeze().sum(axis=0).max()
    else:
        mdlpeak = np.abs(mdl).max()
    if mdlpeak > 0:
        mdl *= amp/mdlpeak

    return mdl


class formaldehyde_mm_model(model.SpectralModel):
    pass

formaldehyde_mm_fitter = formaldehyde_mm_model(formaldehyde_mm, 3,
        parnames=['amp','center','width'],
        parlimited=[(False,False),(False,False), (True,False)],
        parlimits=[(0,0), (0,0), (0,0)],
        shortvarnames=("A","v","\\sigma"), # specify the parameter names (TeX is OK)
        fitunit='Hz' )

formaldehyde_mm_vheight_fitter = formaldehyde_mm_model(fitter.vheightmodel(formaldehyde_mm), 4,
        parnames=['height','amp','center','width'],
        parlimited=[(False,False),(False,False),(False,False), (True,False)],
        parlimits=[(0,0), (0,0), (0,0), (0,0)],
        shortvarnames=("H","A","v","\\sigma"), # specify the parameter names (TeX is OK)
        fitunit='Hz' )


try:
    import pymodelfit

    class pmfFormaldehydeModel(pymodelfit.FunctionModel1DAuto):
        def f(self, x, amp0=1.0, xoff_v0=0.0,width0=1.0):
            return formaldehyde(x,
                    amp=amp0,
                    xoff_v=xoff_v0,width=width0)

    class pmfFormaldehydeModelVtau(pymodelfit.FunctionModel1DAuto):
        def f(self, x, Tex0=1.0, tau0=0.01, xoff_v0=0.0, width0=1.0):
            return formaldehyde_vtau(x,
                    Tex=Tex0, tau=tau0,
                    xoff_v=xoff_v0,width=width0)
except ImportError:
    pass
