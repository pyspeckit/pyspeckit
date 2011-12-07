"""
Fit a line based on parameters output from a grid of RADEX models
"""
import numpy as np
from mpfit import mpfit
from .. import units
from . import fitter,model
import matplotlib.cbook as mpcb
import copy

class radex_model(object):
    def __init__(self, xarr,  
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
            modelfunc=None,
            **kwargs):
        """
        Use a grid of RADEX-computed models to make a model line spectrum

        The RADEX models have to be available somewhere.
        OR they can be passed as arrays.  If as arrays, the form should be:
        texgrid = ((minfreq1,maxfreq1,texgrid1),(minfreq2,maxfreq2,texgrid2))

        xarr must be a SpectroscopicAxis instance
        xoff_v, width are both in km/s.  With is 'sigma'

        grid_vwidth is the velocity assumed when computing the grid in km/s
            this is important because tau = modeltau / width (see, e.g., 
            Draine 2011 textbook pgs 219-230)
        grid_vwidth_scale is True or False: False for LVG, True for Sphere


        A modelfunc must be specified.  Model functions should take an xarr and
        a series of keyword arguments corresponding to the line parameters
        (Tex, tau, xoff_v, and width (gaussian sigma, not FWHM))
        """

        self.modelfunc = modelfunc
        if self.modelfunc is None:
            raise ValueError("Must specify a spectral model function.  See class help for form.")

        if texgrid is None and taugrid is None:
            if path_to_texgrid == '' or path_to_taugrid=='':
                raise IOError("Must specify model grids to use.")
            else:
                self.taugrid = [pyfits.getdata(path_to_taugrid)]
                self.texgrid = [pyfits.getdata(path_to_texgrid)]
                hdr = pyfits.getheader(path_to_taugrid)
                self.yinds,self.xinds = np.indices(self.taugrid[0].shape[1:])
                self.densityarr = (xinds+hdr['CRPIX1']-1)*hdr['CD1_1']+hdr['CRVAL1'] # log density
                self.columnarr  = (yinds+hdr['CRPIX2']-1)*hdr['CD2_2']+hdr['CRVAL2'] # log column
                self.minfreq = (4.8,)
                self.maxfreq = (5.0,)
        elif len(taugrid)==len(texgrid) and hdr is not None:
            self.minfreq,self.maxfreq,self.texgrid = zip(*texgrid)
            self.minfreq,self.maxfreq,self.taugrid = zip(*taugrid)
            self.yinds,self.xinds = np.indices(self.taugrid[0].shape[1:])
            self.densityarr = (xinds+hdr['CRPIX1']-1)*hdr['CD1_1']+hdr['CRVAL1'] # log density
            self.columnarr  = (yinds+hdr['CRPIX2']-1)*hdr['CD2_2']+hdr['CRVAL2'] # log column
        else:
            raise Exception
        
        # Convert X-units to frequency in GHz
        self.xarr = copy.copy(xarr)
        self.xarr.convert_to_unit('Hz', quiet=True)

        #tau = modelgrid.line_params_2D(gridval1,gridval2,densityarr,columnarr,taugrid[temperature_gridnumber,:,:])
        #tex = modelgrid.line_params_2D(gridval1,gridval2,densityarr,columnarr,texgrid[temperature_gridnumber,:,:])

        if debug:
            import pdb; pdb.set_trace()

    def __call__(self, density=4, column=13, xoff_v=0.0, width=1.0,):
        self.gridval1 = np.interp(density, self.densityarr[0,:], xinds[0,:])
        self.gridval2 = np.interp(column, self.columnarr[:,0], yinds[:,0])
        if np.isnan(gridval1) or np.isnan(gridval2):
            raise ValueError("Invalid column/density")

        if scipyOK:
            tau = [scipy.ndimage.map_coordinates(tg[temperature_gridnumber,:,:],np.array([[self.gridval2],[self.gridval1]]),order=1) for tg in self.taugrid]
            tex = [scipy.ndimage.map_coordinates(tg[temperature_gridnumber,:,:],np.array([[self.gridval2],[self.gridval1]]),order=1) for tg in self.texgrid]
        else:
            raise ImportError("Couldn't import scipy, therefore cannot interpolate")

        if verbose:
            print "density %20.12g column %20.12g: tau %20.12g tex %20.12g" % (density, column, tau, tex)
        if debug:
            import pdb; pdb.set_trace()

        return self.modelfunc(self.xarr,Tex=self.tex,tau=tau,xoff_v=xoff_v,width=width,**kwargs)

