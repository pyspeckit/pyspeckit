"""
====================
HCN Hyperfine Fitter
====================
This is an HCN fitter...
ref for line params: http://www.strw.leidenuniv.nl/~moldata/datafiles/hcn@hfs.dat
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
'10-01':88.6339360e9,
'11-01':88.6304160e9,
'12-01':88.6318470e9,
}
aval_dict = {
'10-01':2.4075e-5,
'11-01':2.4075e-5,
'12-01':2.4075e-5,
}
"""
Line strengths of the 15 hyperfine components in J = 1 -  0  transition. The
thickness of the lines indicates their relative weight compared to the others.
Line strengths are normalized in such a way that summing over all initial J = 1
levels gives the degeneracy of the  J = 0 levels, i.e.,  for JF1F = 012,
three for JF1F = 011, and one for JF1F = 010. Thus, the sum over all 15
transitions gives the total spin degeneracy
"""
line_strength_dict = { # effectively the degeneracy per rotation state...
'10-01':1,
'11-01':3,
'12-01':5,
}

line_names = freq_dict.keys()

ckms = units.speedoflight_ms / 1e3 #2.99792458e5
voff_lines_dict = dict([(k,(v-83.6318470e9)/83.6318470e9*ckms) for k,v in freq_dict.iteritems()])

hcn_vtau = hyperfine.hyperfinemodel(line_names, voff_lines_dict, freq_dict, line_strength_dict)
hcn_amp = hcn_vtau.ampfitter
hcn_vtau_fitter = hcn_vtau.fitter
hcn_vtau_vheight_fitter = hcn_vtau.vheight_fitter

def hcn_radex(xarr, density=4, column=13, xoff_v=0.0, width=1.0, 
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

    return hcn_vtau(xarr,Tex=tex,tau=tau,xoff_v=xoff_v,width=width,**kwargs)

