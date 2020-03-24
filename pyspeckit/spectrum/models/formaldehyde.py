"""
===========================
Formaldehyde cm-line fitter
===========================

This is a formaldehyde 1_11-1_10 / 2_12-2_11 fitter.  It includes hyperfine
components of the formaldehyde lines and has both LTE and RADEX LVG based
models

Module API
^^^^^^^^^^
"""
from __future__ import print_function
import numpy as np
from ...mpfit import mpfit
from .. import units
from . import fitter,model,modelgrid
import matplotlib.cbook as mpcb
import copy
from . import hyperfine
from ...specwarnings import warn
from six.moves import xrange
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

try:
    from despotic import cloud
    Democracy=False
except ImportError:
    Democracy=True  # Because it's not despotic :D

from astropy.utils.console import ProgressBar
from astropy.table import Table
import warnings

    
line_names = ['oneone','twotwo','threethree']
line_names = ['oneone_f10', 'oneone_f01', 'oneone_f22', 'oneone_f21',
              'oneone_f12', 'oneone_f11', 'twotwo_f11', 'twotwo_f12',
              'twotwo_f21', 'twotwo_f32', 'twotwo_f33', 'twotwo_f22',
              'twotwo_f23']

# http://adsabs.harvard.edu/abs/1971ApJ...169..429T has the most accurate freqs
# http://adsabs.harvard.edu/abs/1972ApJ...174..463T [twotwo]
central_freq_dict = {
    'oneone':     4.82965996e9,
    'twotwo':     14.48847881e9,
    'threethree': 28.97480e9,
    }
line_strength_dict={
        'oneone_f10':  4.,
        'oneone_f01':  4.,
        'oneone_f22': 15.,
        'oneone_f21':  5.,
        'oneone_f12':  5.,
        'oneone_f11':  3.,
        'twotwo_f11': 15.,
        'twotwo_f12':  5.,
        'twotwo_f21':  5.,
        'twotwo_f32': 5.19,
        'twotwo_f33': 41.49,
        'twotwo_f22': 23.15,
        'twotwo_f23': 5.19,
        'threethree_f22':1,
        'threethree_f44':1,
        'threethree_f33':1,
        }
relative_strength_total_degeneracy={
        'oneone_f10': 36.,
        'oneone_f01': 36.,
        'oneone_f22': 36.,
        'oneone_f21': 36.,
        'oneone_f12': 36.,
        'oneone_f11': 36.,
        'twotwo_f11': 100.01,
        'twotwo_f12': 100.01,
        'twotwo_f21': 100.01,
        'twotwo_f32': 100.01,
        'twotwo_f33': 100.01,
        'twotwo_f22': 100.01,
        'twotwo_f23': 100.01,
        'threethree_f22':3.0,
        'threethree_f44':3.0,
        'threethree_f33':3.0,
        }
hf_freq_dict={
    'oneone_f10':4.82965996e9 - 18.53e3,
    'oneone_f01':4.82965996e9 - 1.34e3,
    'oneone_f22':4.82965996e9 - 0.35e3,
    'oneone_f21':4.82965996e9 + 4.05e3,
    'oneone_f12':4.82965996e9 + 6.48e3,
    'oneone_f11':4.82965996e9 + 11.08e3,
    # 'twotwo_f11':14.48847881e9 - 19.97e3,
    # 'twotwo_f12':14.48847881e9 -  7.03e3,
    # 'twotwo_f21':14.48847881e9 -  2.20e3,
    # 'twotwo_f32':14.48847881e9 +  0.12e3,
    # 'twotwo_f33':14.48847881e9 +  0.89e3,
    # 'twotwo_f22':14.48847881e9 + 10.74e3,
    # 'twotwo_f23':14.48847881e9 + 11.51e3,
    'twotwo_f11':14.48847881e9 - 20.73e3,
    'twotwo_f12':14.48847881e9 -  8.5e3,
    'twotwo_f21':14.48847881e9 -  0.71e3,
    'twotwo_f32':14.48847881e9 +  0.71e3,
    'twotwo_f33':14.48847881e9 +  1.42e3,
    'twotwo_f22':14.48847881e9 + 9.76e3,
    'twotwo_f23':14.48847881e9 + 10.12e3,
    'threethree_f22':28.97478e9,
    'threethree_f44':28.97480e9,
    'threethree_f33':28.97481e9,
}
freq_dict = copy.copy(hf_freq_dict)
freq_dict.update(central_freq_dict)
aval_dict = {
    'oneone':     10**-8.44801,  #64*!pi**4/(3*h*c**3)*nu11**3*mu0**2*(1/2.)
    'twotwo':     10**-7.49373,  #64*!pi**4/(3*h*c**3)*nu22**3*mu0**2*(2/3.)
    'threethree': 10**-6.89179,  #64*!pi**4/(3*h*c**3)*nu33**3*mu0**2*(3/4.)
    }
hf_aval_dict={
        'oneone_f10':10**-8.92509,
        'oneone_f01':10**-8.44797,
        'oneone_f22':10**-8.57294,
        'oneone_f21':10**-9.05004,
        'oneone_f12':10**-8.82819,
        'oneone_f11':10**-9.05009,
        'twotwo_f11':10**-7.61876,
        'twotwo_f12':10**-8.09586,
        'twotwo_f21':10**-8.31771,
        'twotwo_f32':10**-8.44804,
        'twotwo_f33':10**-7.54494,
        'twotwo_f22':10**-7.65221,
        'twotwo_f23':10**-8.30191,
        'threethree_f22':10**-6.94294,
        'threethree_f44':10**-6.91981,
        'threethree_f33':10**-6.96736,
        }
ortho_dict = {
    'oneone':     False,
    'twotwo':     False,
    'threethree': False,
    }
n_ortho = np.arange(0,28,3) # 0..3..27
n_para = np.array([x for x in range(28) if x % 3 != 0])

voff_lines_dict = {
        'oneone': [(hf_freq_dict[f]-freq_dict['oneone'])/freq_dict['oneone']*units.speedoflight_ms for f in hf_freq_dict.keys() if "oneone" in f],
        'twotwo': [(hf_freq_dict[f]-freq_dict['twotwo'])/freq_dict['twotwo']*units.speedoflight_ms for f in hf_freq_dict.keys() if "twotwo" in f],
        'threethree': [(hf_freq_dict[f]-freq_dict['threethree'])/freq_dict['threethree']*units.speedoflight_ms for f in hf_freq_dict.keys() if "threethree" in f],
        }
voff_lines_dict={ # opposite signs of freq offset
        'oneone_f10': + 18.53e3/freq_dict['oneone'] * units.speedoflight_ms / 1000.0,
        'oneone_f01': + 1.34e3 /freq_dict['oneone'] * units.speedoflight_ms / 1000.0,
        'oneone_f22': + 0.35e3 /freq_dict['oneone'] * units.speedoflight_ms / 1000.0,
        'oneone_f21': - 4.05e3 /freq_dict['oneone'] * units.speedoflight_ms / 1000.0,
        'oneone_f12': - 6.48e3 /freq_dict['oneone'] * units.speedoflight_ms / 1000.0,
        'oneone_f11': - 11.08e3/freq_dict['oneone'] * units.speedoflight_ms / 1000.0,
        'twotwo_f11': + 19.97e3/freq_dict['twotwo'] * units.speedoflight_ms / 1000.0,
        'twotwo_f12': +  7.03e3/freq_dict['twotwo'] * units.speedoflight_ms / 1000.0,
        'twotwo_f21': +  2.20e3/freq_dict['twotwo'] * units.speedoflight_ms / 1000.0,
        'twotwo_f32': -  0.12e3/freq_dict['twotwo'] * units.speedoflight_ms / 1000.0,
        'twotwo_f33': -  0.89e3/freq_dict['twotwo'] * units.speedoflight_ms / 1000.0,
        'twotwo_f22': - 10.74e3/freq_dict['twotwo'] * units.speedoflight_ms / 1000.0,
        'twotwo_f23': - 11.51e3/freq_dict['twotwo'] * units.speedoflight_ms / 1000.0,
        'threethree_f22':28.97478e9,
        'threethree_f44':28.97480e9,
        'threethree_f33':28.97481e9,
        }


formaldehyde_vtau = hyperfine.hyperfinemodel(line_names, voff_lines_dict,
                                             freq_dict, line_strength_dict,
                                             relative_strength_total_degeneracy)
formaldehyde_vtau_fitter = formaldehyde_vtau.fitter
formaldehyde_vtau_vheight_fitter = formaldehyde_vtau.vheight_fitter
formaldehyde_vtau_tbg_fitter = formaldehyde_vtau.background_fitter


def build_despotic_grids(gridfile='oh2co_grid_despotic.fits', oh2coAbund=1e-8,
                         nDens=21, logDensLower=2.0, logDensUpper=6.0,
                         nCol=21, logColLower=11.0, logColUpper=15.0,
                         nTemp=51, Tlower=10.0, Tupper=300.0,
                         nDv=5, DvLower=1.0, DvUpper=5.0):
    """
    Generates grids of o-H2CO line intensities using Despotic.  Outputs a astropy Table.
    
    Parameters
    ----------
    gridfile : string
        Name of grid file to output.
    ph2coAbund : float
        Fractional abundance of o-H2CO
    nDens : int
        Number of grid points in the volume density
    logDensLower : float
        log of volume density at lower bound of grid (log(n/cm**-3))
    logDensUpper : float
        log of volume density at upper bound of grid (log(n/cm**-3))
    nCol : int
        Number of grid points in the column density
    logColLower : float
        log of column density of p-H2CO at lower bound of grid (log(N/cm**-2))
    logColUpper : float
        log of column density of p-H2CO at upper bound of grid (log(N/cm**-2))
    nTemp : int
        Number of grid points in the temperature grid
    Tower : float
        temperature at lower bound of grid (K)
    Tupper : float
        temperature at upper bound of grid (K)
    nDv : int
        Number of grid points in the line width
    DvLower : float
        line width (non-thermal) at lower bound of grid (km/s)
    DvUpper : float
        line width (non-thermal) at upper bound of grid (km/s)

    """

    if Democracy:
        raise Exception("No despotic install found.  Cannot build grids")

    
    core = cloud(fileName="protostellarCore.desp", verbose=True)

    nlower = logDensLower
    nupper = logDensUpper

    Nlower = logColLower
    Nupper = logColUpper

    Temps = np.linspace(Tlower, Tupper, nTemp)
    Cols = 1e1**np.linspace(Nlower, Nupper, nCol)
    Densities = 1e1**(np.linspace(nlower, nupper, nDens))
    LineWidth = np.linspace(DvLower, DvUpper, nDv)

    outtable = Table(names = ['Tex_110_111', 'Tex_211_212', 'Tex_312_313',
                              'tau_110_111', 'tau_211_212', 'tau_312_313',
                              'Temperature', 'Column', 'nH2', 'sigmaNT'])

    TempArr, ColArr, DensArr, DvArr = np.meshgrid(Temps,
                                                  Cols,
                                                  Densities,
                                                  LineWidth)

    for T, N, n, dv in ProgressBar(zip(TempArr.flatten(),
                                   ColArr.flatten(),
                                   DensArr.flatten(),
                                   DvArr.flatten())):
        core.colDen = N/oh2coAbund
        core.Tg = T
        core.Td = T
        core.nH = n
        core.sigmaNT = dv
        lines = core.lineLum('o-h2co')
        outtable.add_row()
        outtable[-1]['Tex_110_111'] = lines[0]['Tex']
        outtable[-1]['tau_110_111'] = lines[0]['tau']
        outtable[-1]['Tex_211_212'] = lines[2]['Tex']
        outtable[-1]['tau_211_212'] = lines[2]['tau']
        outtable[-1]['Tex_312_313'] = lines[5]['Tex']
        outtable[-1]['tau_312_313'] = lines[5]['tau']
        outtable[-1]['Temperature'] = T
        outtable[-1]['Column'] = N
        outtable[-1]['nH2'] = n
        outtable[-1]['sigmaNT'] = dv

    outtable.write(gridfile, format='fits', overwrite=True)

    
def formaldehyde_despotic_functions(gridtable):
    """
    This builds interpolation functions for use in fitting.

    Parameters
    ----------
    gridtable : str
       Name of grid in astropy table

    Returns
    -------
    h2co_110_111, h2co_211_212, h2co_312_313 : function
       Functions that return the excitation temperature and optical depth given input density,
       temperature, column density and line width.

    """
    if gridtable is None:
        warnings.warn("No gridfile found.  Building grids using despotic")
        try:
            build_despotic_grids('ph2co_grid_despotic.fits')
            gridtable = Table.read('ph2co_grid_despotic.fits')
        except:  # TODO -- make this more specific
            warnings.warn("Failed to build functions because no grids available")
            return

    DensArr = np.sort(np.unique(gridtable['nH2']))
    ColArr = np.sort(np.unique(gridtable['Column']))
    TempArr = np.sort(np.unique(gridtable['Temperature']))
    DvArr = np.sort(np.unique(gridtable['sigmaNT']))
    GridData_Tex_110_111 = np.zeros((len(DensArr), len(ColArr),
                                     len(TempArr), len(DvArr))) + np.nan
    GridData_Tex_211_212 = np.zeros((len(DensArr), len(ColArr),
                                     len(TempArr), len(DvArr))) + np.nan
    GridData_Tex_312_313 = np.zeros((len(DensArr), len(ColArr),
                                     len(TempArr), len(DvArr))) + np.nan    
    GridData_tau_110_111 = np.zeros((len(DensArr), len(ColArr),
                                     len(TempArr), len(DvArr))) + np.nan
    GridData_tau_211_212 = np.zeros((len(DensArr), len(ColArr),
                                     len(TempArr), len(DvArr))) + np.nan
    GridData_tau_312_313 = np.zeros((len(DensArr), len(ColArr),
                                     len(TempArr), len(DvArr))) + np.nan    

    ii = np.interp(gridtable['nH2'], DensArr, np.arange(len(DensArr))).astype(np.int)
    jj = np.interp(gridtable['Column'], ColArr, np.arange(len(ColArr))).astype(np.int)
    kk = np.interp(gridtable['Temperature'], TempArr, np.arange(len(TempArr))).astype(np.int)
    ll = np.interp(gridtable['sigmaNT'], DvArr, np.arange(len(DvArr))).astype(np.int)

    GridData_Tex_110_111[ii, jj, kk, ll] = gridtable['Tex_110_111']
    GridData_Tex_211_212[ii, jj, kk, ll] = gridtable['Tex_211_212']
    GridData_Tex_312_313[ii, jj, kk, ll] = gridtable['Tex_312_313']
    GridData_tau_110_111[ii, jj, kk, ll] = gridtable['tau_110_111']
    GridData_tau_211_212[ii, jj, kk, ll] = gridtable['tau_211_212']
    GridData_tau_312_313[ii, jj, kk, ll] = gridtable['tau_312_313']

    def h2co_110_111(logdensity=4, logcolumn=13, temperature=25, sigmav=2.0):
        iidx = np.interp(logdensity, np.log10(DensArr), np.arange(len(DensArr)))
        jidx = np.interp(logcolumn, np.log10(ColArr), np.arange(len(ColArr)))
        kidx = np.interp(temperature, TempArr, np.arange(len(TempArr)))
        lidx = np.interp(sigmav, DvArr, np.arange(len(DvArr)))
        xvec = np.array([iidx, jidx, kidx, lidx])
        xvec.shape += (1,)
        Tex = scipy.ndimage.interpolation.map_coordinates(GridData_Tex_110_111,
                                                          xvec)
        tau = scipy.ndimage.interpolation.map_coordinates(GridData_tau_110_111,
                                                          xvec)
        return (Tex, tau)

    def h2co_211_212(logdensity=4, logcolumn=13, temperature=25, sigmav=2.0):
        iidx = np.interp(logdensity, np.log10(DensArr), np.arange(len(DensArr)))
        jidx = np.interp(logcolumn, np.log10(ColArr), np.arange(len(ColArr)))
        kidx = np.interp(temperature, TempArr, np.arange(len(TempArr)))
        lidx = np.interp(sigmav, DvArr, np.arange(len(DvArr)))
        xvec = np.array([iidx, jidx, kidx, lidx])
        xvec.shape += (1,)
        Tex = scipy.ndimage.interpolation.map_coordinates(GridData_Tex_211_212, xvec)
        tau = scipy.ndimage.interpolation.map_coordinates(GridData_tau_211_212, xvec)
        return (Tex, tau)

    def h2co_312_313(logdensity=4, logcolumn=13, temperature=25, sigmav=2.0):
        iidx = np.interp(logdensity, np.log10(DensArr), np.arange(len(DensArr)))
        jidx = np.interp(logcolumn, np.log10(ColArr), np.arange(len(ColArr)))
        kidx = np.interp(temperature, TempArr, np.arange(len(TempArr)))
        lidx = np.interp(sigmav, DvArr, np.arange(len(DvArr)))
        xvec = np.array([iidx, jidx, kidx, lidx])
        xvec.shape += (1,)
        Tex = scipy.ndimage.interpolation.map_coordinates(GridData_Tex_312_313, xvec)
        tau = scipy.ndimage.interpolation.map_coordinates(GridData_tau_312_313, xvec)
        return (Tex, tau)
    return (h2co_110_111, h2co_211_212, h2co_312_313)

def formaldehyde_despotic(xarr, 
                          temperature=25,
                          column=13,
                          density=4,
                          xoff_v=0.0,
                          width=1.0, 
                          grid_vwidth=1.0,
                          h2co_110_111=None,
                          h2co_211_212=None,
                          h2co_312_313=None,
                          debug=False,
                          verbose=False,
                          **kwargs):
    """
    Fitter to o-H2CO using despotic grids.  Requires building grids and passing in 
    functions for interpolating the h2co transition optical depth and 
    excitation temperatures. 
    """
    # The grids are built with a specified NON-THERMAL limit
    width_NT = np.sqrt(width**2 - 1.38e-26 * temperature / (30 * 1.66e-24))

    if np.isnan(width_NT):
        width_NT=0.00
        
    Tex110_111, tau110_111 = h2co_110_111(logdensity=density,
                                          logcolumn=column,
                                          temperature=temperature,
                                          sigmav=width_NT)
    
    Tex211_212, tau211_212 = h2co_211_212(logdensity=density,
                                          logcolumn=column,
                                          temperature=temperature,
                                          sigmav=width_NT)

    Tex312_313, tau312_313 = h2co_312_313(logdensity=density,
                                          logcolumn=column,
                                          temperature=temperature,
                                          sigmav=width_NT)

    tex = [Tex110_111, Tex211_212, Tex312_313]
    tau = [tau110_111, tau211_212, tau312_313]
    minfreq = [3.0, 13.00, 26.00]
    maxfreq = [6.0, 16.00, 30.00]
    spec = np.sum([(formaldehyde_vtau(xarr,Tex=float(tex[ii]),
                                      tau=float(tau[ii]),
                                      xoff_v=xoff_v,
                                      width=width, **kwargs)
                    * (xarr.as_unit('GHz').value > minfreq[ii]) *
                    (xarr.as_unit('GHz').value < maxfreq[ii]))
                   for ii in xrange(len(tex))], 
                  axis=0)
    return spec

def formaldehyde_despotic_twocomp(xarr, 
                                  temperature=25,
                                  column=13,
                                  density=4,
                                  xoff_v=0.0,
                                  width=1.0, 
                                  temperature2=25,
                                  column2=13,
                                  density2=4,
                                  xoff_v2=0.0,
                                  width2=1.0, 
                                  grid_vwidth=1.0,
                                  h2co_110_111=None,
                                  h2co_211_212=None,
                                  h2co_312_313=None,
                                  debug=False,
                                  verbose=False,
                                  **kwargs):
    """
    Fitter to o-H2CO using despotic grids.  Requires building grids and passing in 
    functions for interpolating the h2co transition optical depth and 
    excitation temperatures. 
    """
    # The grids are built with a specified NON-THERMAL limit
    width_NT = np.sqrt(width**2 - 1.38e-26 * temperature / (30 * 1.66e-24))
    width_NT2 = np.sqrt(width2**2 - 1.38e-26 * temperature2 / (30 * 1.66e-24))

    if np.isnan(width_NT):
        width_NT=0.0
    if np.isnan(width_NT2):
        width_NT2=0.0
        
    Tex110_111, tau110_111 = h2co_110_111(logdensity=density,
                                          logcolumn=column,
                                          temperature=temperature,
                                          sigmav=width_NT)
    
    Tex211_212, tau211_212 = h2co_211_212(logdensity=density,
                                          logcolumn=column,
                                          temperature=temperature,
                                          sigmav=width_NT)

    Tex312_313, tau312_313 = h2co_312_313(logdensity=density,
                                          logcolumn=column,
                                          temperature=temperature,
                                          sigmav=width_NT)

    tex = [Tex110_111, Tex211_212, Tex312_313]
    tau = [tau110_111, tau211_212, tau312_313]
    minfreq = [3.0, 13.00, 26.00]
    maxfreq = [6.0, 16.00, 30.00]
    spec = np.sum([(formaldehyde_vtau(xarr,Tex=float(tex[ii]),
                                      tau=float(tau[ii]),
                                      xoff_v=xoff_v,
                                      width=width, **kwargs)
                    * (xarr.as_unit('GHz').value > minfreq[ii]) *
                    (xarr.as_unit('GHz').value < maxfreq[ii]))
                   for ii in xrange(len(tex))], 
                  axis=0)

    Tex110_111, tau110_111 = h2co_110_111(logdensity=density2,
                                          logcolumn=column2,
                                          temperature=temperature2,
                                          sigmav=width_NT2)
    
    Tex211_212, tau211_212 = h2co_211_212(logdensity=density2,
                                          logcolumn=column2,
                                          temperature=temperature2,
                                          sigmav=width_NT2)

    Tex312_313, tau312_313 = h2co_312_313(logdensity=density2,
                                          logcolumn=column2,
                                          temperature=temperature2,
                                          sigmav=width_NT2)

    tex = [Tex110_111, Tex211_212, Tex312_313]
    tau = [tau110_111, tau211_212, tau312_313]
    spec2 = np.sum([(formaldehyde_vtau(xarr,Tex=float(tex[ii]),
                                       tau=float(tau[ii]),
                                       xoff_v=xoff_v2,
                                       width=width2, **kwargs)
                     * (xarr.as_unit('GHz').value > minfreq[ii]) *
                     (xarr.as_unit('GHz').value < maxfreq[ii]))
                    for ii in xrange(len(tex))], 
                   axis=0)

    
    specsum = spec+spec2
    print (temperature, temperature2, column, column2, density, density2, xoff_v2, xoff_v, width_NT, width_NT2)
    return specsum



def formaldehyde_radex(xarr, density=4, column=13, xoff_v=0.0, width=1.0,
                       grid_vwidth=1.0, grid_vwidth_scale=False, texgrid=None,
                       taugrid=None, hdr=None, path_to_texgrid='',
                       path_to_taugrid='', temperature_gridnumber=3,
                       debug=False, verbose=False, **kwargs):
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
    xarr = xarr.as_unit('Hz', quiet=True)

    tau_nu_cumul = np.zeros(len(xarr))

    gridval1 = np.interp(density, densityarr[0,:], xinds[0,:])
    gridval2 = np.interp(column, columnarr[:,0], yinds[:,0])
    if np.isnan(gridval1) or np.isnan(gridval2):
        raise ValueError("Invalid column/density")

    if scipyOK:
        slices = [temperature_gridnumber] + [slice(int(np.floor(gv)),int(np.floor(gv))+2) for gv in (gridval2,gridval1)]
        tau = [scipy.ndimage.map_coordinates(tg[slices],np.array([[gridval2%1],[gridval1%1]]),order=1) for tg in taugrid]
        tex = [scipy.ndimage.map_coordinates(tg[slices],np.array([[gridval2%1],[gridval1%1]]),order=1) for tg in texgrid]
    else:
        raise ImportError("Couldn't import scipy, therefore cannot interpolate")
    #tau = modelgrid.line_params_2D(gridval1,gridval2,densityarr,columnarr,taugrid[temperature_gridnumber,:,:])
    #tex = modelgrid.line_params_2D(gridval1,gridval2,densityarr,columnarr,texgrid[temperature_gridnumber,:,:])

    if verbose:
        print("density %20.12g column %20.12g: tau %20.12g tex %20.12g" % (density, column, tau, tex))

    if debug:
        import pdb; pdb.set_trace()

    spec = np.sum([(formaldehyde_vtau(xarr,Tex=float(tex[ii]),tau=float(tau[ii]),xoff_v=xoff_v,width=width, **kwargs)
                * (xarr.as_unit('GHz')>minfreq[ii]) * (xarr.as_unit('GHz')<maxfreq[ii])) for ii in xrange(len(tex))],
                axis=0)

    return spec

def formaldehyde_radex_orthopara_temp(xarr, density=4, column=13,
                                      orthopara=1.0, temperature=15.0,
                                      xoff_v=0.0, width=1.0,
                                      Tbackground1=2.73,
                                      Tbackground2=2.73,
                                      grid_vwidth=1.0,
                                      grid_vwidth_scale=False, texgrid=None,
                                      taugrid=None, hdr=None,
                                      path_to_texgrid='', path_to_taugrid='',
                                      debug=False, verbose=False,
                                      getpars=False, **kwargs):
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
            minfreq = (4.8,)
            maxfreq = (5.0,)
    elif len(taugrid)==len(texgrid) and hdr is not None:
        minfreq,maxfreq,texgrid = zip(*texgrid)
        minfreq,maxfreq,taugrid = zip(*taugrid)
    else:
        raise Exception

    densityarr = (np.arange(taugrid[0].shape[3])+hdr['CRPIX1']-1)*hdr['CD1_1']+hdr['CRVAL1'] # log density
    columnarr  = (np.arange(taugrid[0].shape[2])+hdr['CRPIX2']-1)*hdr['CD2_2']+hdr['CRVAL2'] # log column
    temparr    = (np.arange(taugrid[0].shape[1])+hdr['CRPIX3']-1)*hdr['CDELT3']+hdr['CRVAL3'] # temperature
    oprarr     = (np.arange(taugrid[0].shape[0])+hdr['CRPIX4']-1)*hdr['CDELT4']+hdr['CRVAL4'] # log ortho/para ratio

    gridval1 = np.interp(density,     densityarr,  np.arange(len(densityarr)))
    gridval2 = np.interp(column,      columnarr,   np.arange(len(columnarr)))
    gridval3 = np.interp(temperature, temparr,     np.arange(len(temparr)))
    gridval4 = np.interp(orthopara,   oprarr,      np.arange(len(oprarr)))
    if np.isnan(gridval1) or np.isnan(gridval2):
        raise ValueError("Invalid column/density")

    if scipyOK:
        slices = [slice(int(np.floor(gv)),int(np.floor(gv)+2))
                  for gv in (gridval4,gridval3,gridval2,gridval1)]
        tau = [scipy.ndimage.map_coordinates(tg[slices],
                                             np.array([[gridval4 % 1],
                                                       [gridval3 % 1],
                                                       [gridval2 % 1],
                                                       [gridval1 % 1]]),
                                             order=1, prefilter=False)
               for tg in taugrid]
        tex = [scipy.ndimage.map_coordinates(tg[slices],
                                             np.array([[gridval4 % 1],
                                                       [gridval3 % 1],
                                                       [gridval2 % 1],
                                                       [gridval1 % 1]]),
                                             order=1,prefilter=False)
               for tg in texgrid]
    else:
        raise ImportError("Couldn't import scipy, therefore cannot interpolate")
    #tau = modelgrid.line_params_2D(gridval1,gridval2,densityarr,columnarr,taugrid[temperature_gridnumber,:,:])
    #tex = modelgrid.line_params_2D(gridval1,gridval2,densityarr,columnarr,texgrid[temperature_gridnumber,:,:])

    # there can be different background temperatures at each frequency
    tbg = [Tbackground1,Tbackground2]

    if verbose:
        print("density %20.12g   column: %20.12g   temperature: %20.12g   opr: %20.12g   xoff_v: %20.12g   width: %20.12g" % (density, column, temperature, orthopara, xoff_v, width))
        print("tau: ",tau," tex: ",tex)
        print("minfreq: ",minfreq," maxfreq: ",maxfreq)
        print("tbg: ",tbg)

    if debug > 1:
        import pdb; pdb.set_trace()

    if getpars:
        return tau,tex

    spec = np.sum([(formaldehyde_vtau(xarr.as_unit('Hz', quiet=True),
                                      Tex=float(tex[ii]), tau=float(tau[ii]),
                                      Tbackground=tbg[ii], xoff_v=xoff_v,
                                      width=width, **kwargs)
                    * (xarr.as_unit('GHz')>minfreq[ii])
                    * (xarr.as_unit('GHz')<maxfreq[ii]))
                   for ii in xrange(len(tex))],
                  axis=0)

    return spec


def formaldehyde(xarr, amp=1.0, xoff_v=0.0, width=1.0,
        return_hyperfine_components=False, texscale=0.01, tau=0.01, **kwargs):
    """
    Generate a model Formaldehyde spectrum based on simple gaussian parameters

    the "amplitude" is an essentially arbitrary parameter; we therefore define
    it to be Tex given tau=0.01 when passing to the fitter
    The final spectrum is then rescaled to that value
    """

    mdl = formaldehyde_vtau(xarr, Tex=amp*texscale, tau=tau, xoff_v=xoff_v,
            width=width,
            return_tau=True,
            return_hyperfine_components=return_hyperfine_components, **kwargs)
    if return_hyperfine_components:
        mdlpeak = np.abs(mdl).squeeze().sum(axis=0).max()
    else:
        mdlpeak = np.abs(mdl).max()
    if mdlpeak > 0:
        mdl *= amp/mdlpeak

    return mdl

def formaldehyde_pyradex(xarr, density=4, column=13, temperature=20,
                         xoff_v=0.0, opr=1.0, width=1.0, tbackground=2.73,
                         grid_vwidth=1.0, debug=False, verbose=False,
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
    """

    raise NotImplementedError("Not done yet.")
    import pyradex

    # Convert X-units to frequency in GHz
    xarr = xarr.as_unit('Hz', quiet=True)

    tb_nu_cumul = np.zeros(len(xarr))

    R = pyradex.Radex(molecule='oh2co-h2', column=column,
                      temperature=temperature, density=10**density,
                      tbackground=tbackground,)

    spec = np.sum([(formaldehyde_vtau(xarr,Tex=float(tex[ii]),tau=float(tau[ii]),xoff_v=xoff_v,width=width, **kwargs)
                * (xarr.as_unit('GHz')>minfreq[ii]) * (xarr.as_unit('GHz')<maxfreq[ii])) for ii in xrange(len(tex))],
                axis=0)

    return spec

class formaldehyde_model(model.SpectralModel):
    def formaldehyde_integral(self, modelpars, linename='oneone'):
        """
        Return the integral of the individual components (ignoring height)
        """
        raise NotImplementedError("Not implemented, but the integral is just amplitude * width * sqrt(2*pi)")
        # produced by directly computing the integral of gaussians and formaldehydeians as a function of
        # line width and then fitting that with a broken logarithmic power law
        # The errors are <0.5% for all widths
        formaldehyde_to_gaussian_ratio_coefs = {
                'lt0.1_oneone': np.array([ -5.784020,-40.058798,-111.172706,-154.256411,-106.593122,-28.933119]),
                'gt0.1_oneone': np.array([  0.038548, -0.071162, -0.045710,  0.183828, -0.145429,  0.040039]),
                'lt0.1_twotwo': np.array([  1.156561,  6.638570, 11.782065, -0.429536,-24.860297,-27.902274, -9.510288]),
                'gt0.1_twotwo': np.array([ -0.090646,  0.078204,  0.123181, -0.175590,  0.089506, -0.034687,  0.008676]),
                }


        integ = 0
        if len(modelpars) % 3 == 0:
            for amp,cen,width in np.reshape(modelpars,[len(modelpars)/3,3]):
                gaussint = amp*width*np.sqrt(2.0*np.pi)
                cftype = "gt0.1_"+linename if width > 0.1 else "lt0.1_"+linename
                correction_factor = 10**np.polyval(formaldehyde_to_gaussian_ratio_coefs[cftype], np.log10(width) )
                # debug statement print("Two components of the integral: amp %g, width %g, gaussint %g, correction_factor %g " % (amp,width,gaussint,correction_factor))
                integ += gaussint*correction_factor

        return integ

formaldehyde_fitter = formaldehyde_model(formaldehyde, 3,
        parnames=['amp','center','width'],
        parlimited=[(False,False),(False,False), (True,False)],
        parlimits=[(0,0), (0,0), (0,0)],
        shortvarnames=("A","v","\\sigma"), # specify the parameter names (TeX is OK)
        fitunit='Hz' )

formaldehyde_vheight_fitter = formaldehyde_model(fitter.vheightmodel(formaldehyde), 4,
        parnames=['height','amp','center','width'],
        parlimited=[(False,False),(False,False),(False,False), (True,False)],
        parlimits=[(0,0), (0,0), (0,0), (0,0)],
        shortvarnames=("H","A","v","\\sigma"), # specify the parameter names (TeX is OK)
        fitunit='Hz' )

# Create a tau-only fit:
def formaldehyde_radex_tau(xarr, density=4, column=13, xoff_v=0.0, width=1.0,
                           grid_vwidth=1.0, grid_vwidth_scale=False,
                           taugrid=None, hdr=None, path_to_taugrid='',
                           temperature_gridnumber=3, debug=False,
                           verbose=False, return_hyperfine_components=False,
                           **kwargs):
    """
    Use a grid of RADEX-computed models to make a model line spectrum

     * uses hyperfine components
     * assumes *tau* varies but *tex* does not!

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

    if verbose:
        print("Parameters: dens=%f, column=%f, xoff=%f, width=%f" % (density, column, xoff_v, width))

    if taugrid is None:
        if path_to_taugrid=='':
            raise IOError("Must specify model grids to use.")
        else:
            taugrid = [pyfits.getdata(path_to_taugrid)]
            hdr = pyfits.getheader(path_to_taugrid)
            yinds,xinds = np.indices(taugrid[0].shape[1:])
            densityarr = (xinds+hdr['CRPIX1']-1)*hdr['CD1_1']+hdr['CRVAL1'] # log density
            columnarr  = (yinds+hdr['CRPIX2']-1)*hdr['CD2_2']+hdr['CRVAL2'] # log column
            minfreq = (4.8,)
            maxfreq = (5.0,)
    elif hdr is not None:
        minfreq,maxfreq,taugrid = zip(*taugrid)
        yinds,xinds = np.indices(taugrid[0].shape[1:])
        densityarr = (xinds+hdr['CRPIX1']-1)*hdr['CD1_1']+hdr['CRVAL1'] # log density
        columnarr  = (yinds+hdr['CRPIX2']-1)*hdr['CD2_2']+hdr['CRVAL2'] # log column
    else:
        raise Exception

    # Convert X-units to frequency in GHz
    xarr = xarr.as_unit('Hz', quiet=True)

    gridval1 = np.interp(density, densityarr[0,:], xinds[0,:])
    gridval2 = np.interp(column, columnarr[:,0], yinds[:,0])
    if np.isnan(gridval1) or np.isnan(gridval2):
        raise ValueError("Invalid column/density")

    if scipyOK:
        slices = [temperature_gridnumber] + [slice(np.floor(gv),np.floor(gv)+2) for gv in (gridval2,gridval1)]
        tau = [scipy.ndimage.map_coordinates(tg[slices],np.array([[gridval2%1],[gridval1%1]]),order=1) for tg in taugrid]
    else:
        raise ImportError("Couldn't import scipy, therefore cannot interpolate")

    # let the hyperfine module determine the hyperfine components, and pass all of them here
    spec_components = [(formaldehyde_vtau(xarr.as_unit('Hz', quiet=True),
        tau=float(tau[ii]), xoff_v=xoff_v, width=width,
        return_tau=True, return_hyperfine_components=True, **kwargs) *
            (xarr.as_unit('GHz')>minfreq[ii]) *
            (xarr.as_unit('GHz')<maxfreq[ii]))
                for ii in  xrange(len(tau))]

    # get an array of [n_lines, n_hyperfine, len(xarr)]
    if return_hyperfine_components:
        return np.array(spec_components).sum(axis=0)
    else:
        return np.sum(spec_components, axis=0).sum(axis=0)


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
