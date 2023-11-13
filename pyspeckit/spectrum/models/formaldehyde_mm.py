"""
===========================
Formaldehyde mm-line fitter
===========================

This is a formaldehyde 3_03-2_02 / 3_22-221 and 3_03-2_02/3_21-2_20 fitter.
It is based entirely on RADEX models.

Module API
^^^^^^^^^^
"""
from __future__ import print_function
import numpy as np
from . import hyperfine
from . import fitter,model#,modelgrid
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

def build_despotic_grids(gridfile='ph2co_grid_despotic.fits', ph2coAbund=1e-8,
                         nDens=21, logDensLower=2.0, logDensUpper=6.0,
                         nCol=21, logColLower=11.0, logColUpper=15.0,
                         nTemp=51, Tlower=10.0, Tupper=300.0,
                         nDv=5, DvLower=1.0, DvUpper=5.0):
    """
    Generates grids of p-H2CO line intensities using Despotic.  Outputs a astropy Table.

    Parameters
    ----------
    gridfile : string
        Name of grid file to output.
    ph2coAbund : float
        Fractional abundance of p-H2CO
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

    outtable = Table(names = ['Tex_303_202', 'Tex_322_221', 'Tex_321_220',
                              'tau_303_202', 'tau_322_221', 'tau_321_220',
                              'Temperature', 'Column', 'nH2', 'sigmaNT'])

    TempArr, ColArr, DensArr, DvArr = np.meshgrid(Temps,
                                                  Cols,
                                                  Densities,
                                                  LineWidth)

    for T, N, n, dv in ProgressBar(zip(TempArr.flatten(),
                                   ColArr.flatten(),
                                   DensArr.flatten(),
                                   DvArr.flatten())):
        core.colDen = N/ph2coAbund
        core.Tg = T
        core.Td = T
        core.nH = n
        core.sigmaNT = dv
        lines = core.lineLum('p-h2co')
        outtable.add_row()
        outtable[-1]['Tex_303_202'] = lines[2]['Tex']
        outtable[-1]['tau_303_202'] = lines[2]['tau']
        outtable[-1]['Tex_322_221'] = lines[9]['Tex']
        outtable[-1]['tau_322_221'] = lines[9]['tau']
        outtable[-1]['Tex_321_220'] = lines[12]['Tex']
        outtable[-1]['tau_321_220'] = lines[12]['tau']
        outtable[-1]['Temperature'] = T
        outtable[-1]['Column'] = N
        outtable[-1]['nH2'] = n
        outtable[-1]['sigmaNT'] = dv

    outtable.write(gridfile, format='fits',overwrite=True)



def formaldehyde_mm_despotic_functions(gridtable):
    """
    This builds interpolation functions for use in fitting.

    Parameters
    ----------
    gridtable : str
       Name of grid in astropy table

    Returns
    -------
    h2co_303_202, h2co_322_221, h2co_321_220 : function
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
    GridData_Tex_303_202 = np.zeros((len(DensArr), len(ColArr),
                                     len(TempArr), len(DvArr))) + np.nan
    GridData_Tex_322_221 = np.zeros((len(DensArr), len(ColArr),
                                     len(TempArr), len(DvArr))) + np.nan
    GridData_Tex_321_220 = np.zeros((len(DensArr), len(ColArr),
                                     len(TempArr), len(DvArr))) + np.nan
    GridData_tau_303_202 = np.zeros((len(DensArr), len(ColArr),
                                     len(TempArr), len(DvArr))) + np.nan
    GridData_tau_322_221 = np.zeros((len(DensArr), len(ColArr),
                                     len(TempArr), len(DvArr))) + np.nan
    GridData_tau_321_220 = np.zeros((len(DensArr), len(ColArr),
                                     len(TempArr), len(DvArr))) + np.nan

    ii = np.interp(gridtable['nH2'], DensArr, np.arange(len(DensArr))).astype(int)
    jj = np.interp(gridtable['Column'], ColArr, np.arange(len(ColArr))).astype(int)
    kk = np.interp(gridtable['Temperature'], TempArr, np.arange(len(TempArr))).astype(int)
    ll = np.interp(gridtable['sigmaNT'], DvArr, np.arange(len(DvArr))).astype(int)

    GridData_Tex_303_202[ii, jj, kk, ll] = gridtable['Tex_303_202']
    GridData_Tex_322_221[ii, jj, kk, ll] = gridtable['Tex_322_221']
    GridData_Tex_321_220[ii, jj, kk, ll] = gridtable['Tex_321_220']
    GridData_tau_303_202[ii, jj, kk, ll] = gridtable['tau_303_202']
    GridData_tau_322_221[ii, jj, kk, ll] = gridtable['tau_322_221']
    GridData_tau_321_220[ii, jj, kk, ll] = gridtable['tau_321_220']

    def h2co_303_202(logdensity=4, logcolumn=13, temperature=25, sigmav=2.0):
        iidx = np.interp(logdensity, np.log10(DensArr), np.arange(len(DensArr)))
        jidx = np.interp(logcolumn, np.log10(ColArr), np.arange(len(ColArr)))
        kidx = np.interp(temperature, TempArr, np.arange(len(TempArr)))
        lidx = np.interp(sigmav, DvArr, np.arange(len(DvArr)))
        xvec = np.array([iidx, jidx, kidx, lidx])
        xvec.shape += (1,)
        Tex = scipy.ndimage.interpolation.map_coordinates(GridData_Tex_303_202,
                                                          xvec)
        tau = scipy.ndimage.interpolation.map_coordinates(GridData_tau_303_202,
                                                          xvec)
        return (Tex, tau)

    def h2co_322_221(logdensity=4, logcolumn=13, temperature=25, sigmav=2.0):
        iidx = np.interp(logdensity, np.log10(DensArr), np.arange(len(DensArr)))
        jidx = np.interp(logcolumn, np.log10(ColArr), np.arange(len(ColArr)))
        kidx = np.interp(temperature, TempArr, np.arange(len(TempArr)))
        lidx = np.interp(sigmav, DvArr, np.arange(len(DvArr)))
        xvec = np.array([iidx, jidx, kidx, lidx])
        xvec.shape += (1,)
        Tex = scipy.ndimage.interpolation.map_coordinates(GridData_Tex_322_221, xvec)
        tau = scipy.ndimage.interpolation.map_coordinates(GridData_tau_322_221, xvec)
        return (Tex, tau)

    def h2co_321_220(logdensity=4, logcolumn=13, temperature=25, sigmav=2.0):
        iidx = np.interp(logdensity, np.log10(DensArr), np.arange(len(DensArr)))
        jidx = np.interp(logcolumn, np.log10(ColArr), np.arange(len(ColArr)))
        kidx = np.interp(temperature, TempArr, np.arange(len(TempArr)))
        lidx = np.interp(sigmav, DvArr, np.arange(len(DvArr)))
        xvec = np.array([iidx, jidx, kidx, lidx])
        xvec.shape += (1,)
        Tex = scipy.ndimage.interpolation.map_coordinates(GridData_Tex_321_220, xvec)
        tau = scipy.ndimage.interpolation.map_coordinates(GridData_tau_321_220, xvec)
        return (Tex, tau)
    return (h2co_303_202, h2co_322_221, h2co_321_220)

def formaldehyde_mm_despotic(xarr,
                             temperature=25,
                             column=13,
                             density=4,
                             xoff_v=0.0,
                             width=1.0,
                             grid_vwidth=1.0,
                             h2co_303_202=None,
                             h2co_322_221=None,
                             h2co_321_220=None,
                             debug=False,
                             verbose=False,
                             **kwargs):
    """
    Fitter to p-H2CO using despotic grids.  Requires building grids and passing in
    functions for interpolating the h2co transition optical depth and
    excitation temperatures.
    """

    Tex303_202, tau303_202 = h2co_303_202(logdensity=density,
                                          logcolumn=column,
                                          temperature=temperature,
                                          sigmav=width)

    Tex322_221, tau322_221 = h2co_322_221(logdensity=density,
                                          logcolumn=column,
                                          temperature=temperature,
                                          sigmav=width)

    Tex321_220, tau321_220 = h2co_321_220(logdensity=density,
                                          logcolumn=column,
                                          temperature=temperature,
                                          sigmav=width)

    tex = [Tex303_202, Tex322_221, Tex321_220]
    tau = [tau303_202, tau322_221, tau321_220]
    minfreq = [218.15, 218.40, 218.7]
    maxfreq = [218.25, 218.55, 218.8]
    spec = np.sum([
        (formaldehyde_mm_vtau(xarr, Tex=float(tex[ii]), tau=float(tau[ii]),
                              xoff_v=xoff_v, width=width, **kwargs)
         * (xarr.as_unit('GHz').value>minfreq[ii]) *
         (xarr.as_unit('GHz').value<maxfreq[ii])) for ii in xrange(len(tex))],
                  axis=0)
    return spec


def formaldehyde_mm_radex(xarr,
        temperature=25,
        column=13,
        density=4,
        xoff_v=0.0,
        width=1.0,
        grid_vwidth=1.0,
        texgrid=None,
        taugrid=None,
        hdr=None,
        path_to_texgrid='',
        path_to_taugrid='',
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

    if texgrid is None and taugrid is None:
        if path_to_texgrid == '' or path_to_taugrid=='':
            raise IOError("Must specify model grids to use.")
        else:
            taugrid = [pyfits.getdata(path_to_taugrid)]
            texgrid = [pyfits.getdata(path_to_texgrid)]
            hdr = pyfits.getheader(path_to_taugrid)
            zinds,yinds,xinds = np.indices(taugrid[0].shape)
            if 'CD1_1' in hdr:
                cd11 = 'CD1_1'
                cd22 = 'CD2_2'
            else:
                cd11 = 'CDELT1'
                cd22 = 'CDELT2'
            densityarr = (xinds+hdr['CRPIX1']-1)*hdr[cd11]+hdr['CRVAL1'] # log density
            columnarr  = (yinds+hdr['CRPIX2']-1)*hdr[cd22]+hdr['CRVAL2'] # log column
            temparr    = (zinds+hdr['CRPIX3']-1)*hdr['CDELT3']+hdr['CRVAL3'] # lin temperature
            minfreq = (218.,)
            maxfreq = (219.,)
    elif len(taugrid)==len(texgrid) and hdr is not None:
        minfreq,maxfreq,texgrid = zip(*texgrid)
        minfreq,maxfreq,taugrid = zip(*taugrid)
        zinds,yinds,xinds = np.indices(taugrid[0].shape)
        if 'CD1_1' in hdr:
            cd11 = 'CD1_1'
            cd22 = 'CD2_2'
        else:
            cd11 = 'CDELT1'
            cd22 = 'CDELT2'
        densityarr = (xinds+hdr['CRPIX1']-1)*hdr[cd11]+hdr['CRVAL1'] # log density
        columnarr  = (yinds+hdr['CRPIX2']-1)*hdr[cd22]+hdr['CRVAL2'] # log column
        temparr    = (zinds+hdr['CRPIX3']-1)*hdr['CDELT3']+hdr['CRVAL3'] # lin temperature
    else:
        raise Exception

    # Convert X-units to frequency in GHz
    xarr = xarr.as_unit('Hz', quiet=True)

    #tau_nu_cumul = np.zeros(len(xarr))

    gridval1 = np.interp(density, densityarr[0,0,:], xinds[0,0,:])
    gridval2 = np.interp(column, columnarr[0,:,0], yinds[0,:,0])
    gridval3 = np.interp(temperature, temparr[:,0,0], zinds[:,0,0])
    if np.isnan(gridval1) or np.isnan(gridval2) or np.isnan(gridval3):
        raise ValueError("Invalid column/density")

    if scipyOK:
        # this is mostly a trick for speed: slice so you only have two thin layers to interpolate
        # between
        #slices = [density_gridnumber] + [slice(np.floor(gv),np.floor(gv)+2) for gv in (gridval2,gridval1)]
        slices = [slice(np.floor(gridval3),np.floor(gridval3)+2),
                  slice(np.floor(gridval2),np.floor(gridval2)+2),
                  slice(np.floor(gridval1),np.floor(gridval1)+2)
                  ]
        tau = [scipy.ndimage.map_coordinates(tg[slices],
            np.array([[gridval3%1],[gridval2%1],[gridval1%1]]),order=1) for tg in taugrid]
        tex = [scipy.ndimage.map_coordinates(tg[slices],
            np.array([[gridval3%1],[gridval2%1],[gridval1%1]]),order=1) for tg in texgrid]
    else:
        raise ImportError("Couldn't import scipy, therefore cannot interpolate")
    #tau = modelgrid.line_params_2D(gridval1,gridval2,densityarr,columnarr,taugrid[temperature_gridnumber,:,:])
    #tex = modelgrid.line_params_2D(gridval1,gridval2,densityarr,columnarr,texgrid[temperature_gridnumber,:,:])

    if verbose:
        for ta,tk in zip(tau,tex):
            print("density %20.12g temperature %20.12g column %20.12g: tau %20.12g tex %20.12g" % (density, temperature, column, ta, tk))

    if debug:
        import pdb; pdb.set_trace()

    spec = np.sum([
        (formaldehyde_mm_vtau(xarr, Tex=float(tex[ii]), tau=float(tau[ii]),
            xoff_v=xoff_v, width=width, **kwargs)
        * (xarr.as_unit('GHz').value>minfreq[ii]) *
         (xarr.as_unit('GHz').value<maxfreq[ii])) for ii in xrange(len(tex))],
                axis=0)

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
