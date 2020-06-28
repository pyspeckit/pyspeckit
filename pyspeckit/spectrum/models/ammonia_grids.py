from astropy.io import fits
from astropy.table import Table
import scipy.ndimage.interpolation as ndinterp
import scipy.interpolate as interp
import numpy as np
import pdb

default_filename = None
sig2fwhm = np.sqrt(8 * np.log(2))

loglogdict = {'oneone': [-1.81326302e-01,  3.34381710e-01],
              'twotwo': [-2.69968439e-02,  1.11981158e-01],
              'threethree': [-2.68340072e-02,  6.72730044e-02],
              'fourfour': [-7.89713752e-05,  2.91441632e-02],
              'fivefive': [-3.25164539e-05,  1.93889654e-02],
              'sixsix': [-1.76846051e-05,  1.38318044e-02],
              'sevenseven': [-1.18320819e-05,  1.03632255e-02],
              'eighteight': [-9.33491439e-06,  8.05543087e-03],
              'ninenine': [ 0.00000000e+00,  0.00000000e+00]}


def line_scaling_function(sigmav, linename, method='loglogfit'):

    # Scaling of intrinsic line width to Gaussian equivalent line
    # width for use in RADEX
    
    if method == 'none':
        return(sigmav)

    if method == 'sig2fwhm':
        return(sigmav * sig2fwhm)
    
    if method == 'loglogfit':
        coeffs = loglogdict[linename]
        return(sigmav *  1e1**(np.log10(sigmav)
                               * coeffs[0] + coeffs[1]))
    
    
def parbounds(gridfile=None):
    # Return boundaries of grid to bound fit parameters
    
    if gridfile is None:
        return(None)

    t = Table.read(gridfile)
    nOutputs = len(t.keys())-4 # this is the number of outputs from the grid.
    modelkeys = (t.keys())[0:nOutputs]
    density_values = np.unique(t['nH2'].data)
    temperature_values = np.unique(t['Temperature'].data)
    column_values = np.unique(t['Column'].data)
    sigmav_values = np.unique(t['sigmav'].data)
    bound_dict = {'logdens': [np.nanmin(np.log10(density_values)),
                              np.nanmax(np.log10(density_values))],
                  'tkin':[np.nanmin(temperature_values),
                          np.nanmax(temperature_values)],
                  'ntot':[np.nanmin(np.log10(column_values)),
                          np.nanmax(np.log10(column_values))],
                  'sigmav':[np.nanmin(sigmav_values),
                            np.nanmax(sigmav_values)]}
    return(bound_dict)

    
def ammonia_grids(gridfile=None):
    if gridfile is None:
        gridfile = default_filename
    try:
        t = Table.read(gridfile)
    except:
        return None
    
    # Given inputs, fold table into a grid.  Assume uniform spacing in
    # log space.

    nOutputs = len(t.keys())-4 # this is the number of outputs from the grid.
    modelkeys = (t.keys())[0:nOutputs]
    density_values = np.unique(t['nH2'].data)
    temperature_values = np.unique(t['Temperature'].data)
    column_values = np.unique(t['Column'].data)
    sigmav_values = np.unique(t['sigmav'].data)

    # Should this be initialized to NaNs instead?  Could be risky.
    grid = np.zeros((nOutputs,density_values.size,temperature_values.size,
                     column_values.size,sigmav_values.size))
    idx1 = np.digitize(t['nH2'].data,density_values)-1
    idx2 = np.digitize(t['Temperature'].data,temperature_values)-1
    idx3 = np.digitize(t['Column'].data,column_values)-1
    idx4 = np.digitize(t['sigmav'].data,sigmav_values)-1

    f1 = interp.interp1d(np.log10(density_values),
                         np.arange(density_values.size),
                         bounds_error=False)
    f2 = interp.interp1d(np.log10(temperature_values),
                         np.arange(temperature_values.size),
                         bounds_error=False)
    f3 = interp.interp1d(np.log10(column_values),
                         np.arange(column_values.size),
                         bounds_error=False)
    f4 = interp.interp1d(np.log10(sigmav_values),
                         np.arange(sigmav_values.size),
                         bounds_error = False)

    # Pack that grid
    for idx,key in enumerate(modelkeys):
                         grid[idx,idx1,idx2,idx3,idx4] = t[key].data
    
    def ammonia_sampler(linename='oneone',
                        logdens=4,
                        tkin=15,
                        ntot=14,
                        fortho=0.0,
                        sigma=0.3,
                        order=2,
                        scaling_method='loglogfit'):
        
        dV = line_scaling_function(sigma, linename,
                                   method=scaling_method)

        param_interps = (np.array(np.c_[0,f1(logdens+np.log10(1-fortho)),
                                       f2(np.log10(tkin)),
                                       f3(ntot),
                                       f4(np.log10(dV))])
                         * np.ones((nOutputs,1)))

        param_interps[:,0] = np.arange(nOutputs)
        interp_params = ndinterp.map_coordinates(grid, param_interps.T,
                                               order=order,cval=np.nan)
        outputdict = {}
        for idx,key in enumerate(modelkeys):
            outputdict[key] = interp_params[idx]

        return outputdict

    return ammonia_sampler
