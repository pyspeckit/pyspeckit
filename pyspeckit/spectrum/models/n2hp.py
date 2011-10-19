"""
===========
N2H+ fitter
===========
Reference for line params: 
Daniel, F., Dubernet, M.-L., Meuwly, M., Cernicharo, J., Pagani, L. 2005, MNRAS 363, 1083
http://www.strw.leidenuniv.nl/~moldata/N2H+.html
http://adsabs.harvard.edu/abs/2005MNRAS.363.1083D

"""
import numpy as np
from mpfit import mpfit
from .. import units
from . import fitter,model,modelgrid
import matplotlib.cbook as mpcb
import copy
try: # for model grid reading
    import pyfits
    pyfitsOK = True
except:
    pyfitsOK = False
try:
    import scipy.interpolate
    import scipy.ndimage
    scipyOK = True
except ImportError:
    scipyOK=False
import hyperfine


freq_dict={
'110-011':93.171617e9,
'112-011':93.171913e9,
'112-012':93.171913e9,
'111-010':93.172048e9,
'111-011':93.172048e9,
'111-012':93.172048e9,
'122-011':93.173475e9,
'122-012':93.173475e9,
'123-012':93.173772e9,
'121-010':93.173963e9,
'121-011':93.173963e9,
'121-012':93.173963e9,
'101-010':93.176261e9,
'101-011':93.176261e9,
'101-012':93.176261e9,
}
aval_dict = {
'110-011':3.628,
'112-011':0.907,
'112-012':2.721,
'111-010':1.209,
'111-011':0.907,
'111-012':1.512,
'122-011':2.721,
'122-012':0.907,
'123-012':3.628,
'121-010':2.015,
'121-011':1.512,
'121-012':0.101,
'101-010':0.403,
'101-011':1.209,
'101-012':2.016,
}
"""
Line strengths of the 15 hyperfine components in J=1-0  transition. The
thickness of the lines indicates their relative weight compared to the others.
Line strengths are normalized in such a way that summing over all initial J = 1
levels gives the degeneracy of the J = 0 levels, i.e.,  for JF1F  012,
three for JF1F  011, and one for JF1F  010. Thus, the sum over all 15
transitions gives the total spin degeneracy
"""
line_strength_dict = { # effectively the degeneracy per rotation state...
'110-011':0.333,
'112-011':0.417,
'112-012':1.250,
'111-010':0.333,
'111-011':0.250,
'111-012':0.417,
'122-011':1.250,
'122-012':0.417,
'123-012':2.330,
'121-010':0.556,
'121-011':0.417,
'121-012':0.028,
'101-010':0.111,
'101-011':0.333,
'101-012':0.55,
}

line_names = freq_dict.keys()

ckms = units.speedoflight_ms / 1e3 #2.99792458e5
voff_lines_dict = dict([(k,(v-93.176261e9)/93.176261e9*ckms) for k,v in freq_dict.iteritems()])

n2hp_vtau = hyperfine.hyperfinemodel(line_names, voff_lines_dict, freq_dict, line_strength_dict)
n2hp_vtau_fitter = n2hp_vtau.fitter
n2hp_vtau_vheight_fitter = n2hp_vtau.vheight_fitter

def n2hp_radex(xarr, density=4, column=13, xoff_v=0.0, width=1.0, 
        grid_vwidth=1.0,
        grid_vwidth_scale=False,
        texgrid=None,
        taugrid=None,
        hdr=None,
        path_to_texgrid='',
        path_to_taugrid='',
        temperature_gridnumber=3,
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

    grid_vwidth is the velocity assumed when computing the grid in km/s
        this is important because tau = modeltau / width (see, e.g., 
        Draine 2011 textbook pgs 219-230)
    grid_vwidth_scale is True or False: False for LVG, True for Sphere
    """

    if texgrid is None and taugrid is None:
        if path_to_texgrid == '' or path_to_taugrid=='':
            raise IOError("Must specify model grids to use.")
        else:
            taugrid = [pyfits.getdata(path_to_taugrid)]
            texgrid = [pyfits.getdata(path_to_texgrid)]
            hdr = pyfits.getheader(path_to_taugrid)
            yinds,xinds = np.indices(taugrid[0].shape[1:])
            densityarr = (xinds+hdr['CRPIX1']-1)*hdr['CD1_1']+hdr['CRVAL1'] # log density
            columnarr  = (yinds+hdr['CRPIX2']-1)*hdr['CD2_2']+hdr['CRVAL2'] # log column
            minfreq = (4.8,)
            maxfreq = (5.0,)
    elif len(taugrid)==len(texgrid) and hdr is not None:
        minfreq,maxfreq,texgrid = zip(*texgrid)
        minfreq,maxfreq,taugrid = zip(*taugrid)
        yinds,xinds = np.indices(taugrid[0].shape[1:])
        densityarr = (xinds+hdr['CRPIX1']-1)*hdr['CD1_1']+hdr['CRVAL1'] # log density
        columnarr  = (yinds+hdr['CRPIX2']-1)*hdr['CD2_2']+hdr['CRVAL2'] # log column
    else:
        raise Exception
    
    # Convert X-units to frequency in GHz
    xarr = copy.copy(xarr)
    xarr.convert_to_unit('Hz', quiet=True)

    tau_nu_cumul = np.zeros(len(xarr))

    gridval1 = np.interp(density, densityarr[0,:], xinds[0,:])
    gridval2 = np.interp(column, columnarr[:,0], yinds[:,0])
    if np.isnan(gridval1) or np.isnan(gridval2):
        raise ValueError("Invalid column/density")

    if scipyOK:
        tau = [scipy.ndimage.map_coordinates(tg[temperature_gridnumber,:,:],np.array([[gridval2],[gridval1]]),order=1) for tg in taugrid]
        tex = [scipy.ndimage.map_coordinates(tg[temperature_gridnumber,:,:],np.array([[gridval2],[gridval1]]),order=1) for tg in texgrid]
    else:
        raise ImportError("Couldn't import scipy, therefore cannot interpolate")
    #tau = modelgrid.line_params_2D(gridval1,gridval2,densityarr,columnarr,taugrid[temperature_gridnumber,:,:])
    #tex = modelgrid.line_params_2D(gridval1,gridval2,densityarr,columnarr,texgrid[temperature_gridnumber,:,:])

    if verbose:
        print "density %20.12g column %20.12g: tau %20.12g tex %20.12g" % (density, column, tau, tex)

    if debug:
        import pdb; pdb.set_trace()

    return n2hp_vtau(xarr,Tex=tex,tau=tau,xoff_v=xoff_v,width=width,**kwargs)

