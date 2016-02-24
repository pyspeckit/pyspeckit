from astropy.io import fits
from astropy.table import Table
import scipy.ndimage.interpolation as ndinterp
import scipy.interpolate as interp
import numpy as np
import pdb

from astropy.utils.data import get_pkg_data_filename
filename = get_pkg_data_filename('data/nh3grid.fits.gz',
                                 package='pyspeckit')

line_scaling = {'oneone':[1.13140205e-01,-4.41680966e-01,-1.48561577e-04],
                'twotwo':[ 0.07469505,-0.34594756,-0.28113087],
                'threethree':[ 0.07081068,-0.34400585,-0.32480547],
                'fourfour':[ 0.07423372,-0.32940902,-0.37835741]}

def ammonia_grids():
    try:
        t = Table.read(filename)
    except:
        return None
    
    # Given inputs, fold table into a grid.  Assume uniform spacing in
    # log space.

    nOutputs = len(t.keys())-4 # this is the number of outputs from the grid.
    modelkeys = (t.keys())[0:nOutputs]
    density_values = np.unique(t['nH2'].data)
    temperature_values = np.unique(t['Temperature'].data)
    column_values = np.unique(t['Column'].data)
    fwhm_values = np.unique(t['FWHM'].data)

    # Should this be initialized to NaNs instead?  Could be risky.
    grid = np.zeros((nOutputs,density_values.size,temperature_values.size,
                     column_values.size,fwhm_values.size))
    idx1 = np.digitize(t['nH2'].data,density_values)-1
    idx2 = np.digitize(t['Temperature'].data,temperature_values)-1
    idx3 = np.digitize(t['Column'].data,column_values)-1
    idx4 = np.digitize(t['FWHM'].data,fwhm_values)-1

    f1 = interp.interp1d(np.log10(density_values),
                         np.arange(density_values.size),
                         bounds_error=False)
    f2 = interp.interp1d(np.log10(temperature_values),
                         np.arange(temperature_values.size),
                         bounds_error=False)
    f3 = interp.interp1d(np.log10(column_values),
                         np.arange(column_values.size),
                         bounds_error=False)
    f4 = interp.interp1d(np.log10(fwhm_values),
                         np.arange(fwhm_values.size),
                         bounds_error = False)

    # Pack that grid
    for idx,key in enumerate(modelkeys):
                         grid[idx,idx1,idx2,idx3,idx4] = t[key].data
    
    def ammonia_sampler(density = 1e4, temperature = 15,
                        column_density =1e14, fortho = 0.0, sigma = 0.3):

#        linelist = ['oneone','twotwo','threethree','fourfour']
        fwhm = sigma * (8*np.log(2))**0.5 
        #  np.zeros(len(linelist))
#        for idx, linename in enumerate(linelist):
#            fwhm[idx] = 1e1**np.polyval(line_scaling[linename],np.log10(sigma))

        param_interps = np.array(np.c_[0,f1(np.log10(density*(1-fortho))),\
                                        f2(np.log10(temperature)),\
                                        f3(np.log10(column_density)),\
                                        f4(np.log10(fwhm))])*\
                                        np.ones((nOutputs,1))

        param_interps[:,0] = np.arange(nOutputs)
        interp_params = ndinterp.map_coordinates(grid, param_interps.T,
                                               order = 2,cval=np.nan)
        outputdict = {}
        for idx,key in enumerate(modelkeys):
            outputdict[key] = interp_params[idx]

        return outputdict

    return ammonia_sampler
