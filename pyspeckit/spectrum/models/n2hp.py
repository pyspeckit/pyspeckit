"""
===========
N2H+ fitter
===========
Reference for line params: 
Dore (Private Communication), improving on the determinations from 
L. Pagani, F. Daniel, and M. L. Dubernet A&A 494, 719-727 (2009)
DOI: 10.1051/0004-6361:200810570

http://www.strw.leidenuniv.nl/~moldata/N2H+.html

http://adsabs.harvard.edu/abs/2005MNRAS.363.1083D

"""
from __future__ import print_function
import numpy as np
import matplotlib.cbook as mpcb
import copy
try:
    from astropy.io import fits as pyfits
except ImportError:
    import pyfits
try:
    import scipy.interpolate
    import scipy.ndimage
    scipyOK = True
except ImportError:
    scipyOK=False

from ...mpfit import mpfit
from .. import units
from . import fitter,model,modelgrid
from . import hyperfine
import astropy.units as u

freq_dict_cen ={
                'J1-0':  93173.7637e6,
                'J2-1': 186344.8420e6,
                'J3-2': 279511.8325e6,
               }

voff_lines_dict={
    ####### J 1-0
    'J1-0_01': -7.9930,
    'J1-0_02': -7.9930,
    'J1-0_03': -7.9930,
    'J1-0_04': -0.6112,
    'J1-0_05': -0.6112,
    'J1-0_06': -0.6112,
    'J1-0_07': 0.0000,
    'J1-0_08': 0.9533,
    'J1-0_09': 0.9533,
    'J1-0_10': 5.5371,
    'J1-0_11': 5.5371,
    'J1-0_12': 5.5371,
    'J1-0_13': 5.9704,
    'J1-0_14': 5.9704,
    'J1-0_15': 6.9238,
    ####### J 2-1
    'J2-1_01': -4.6258,
    'J2-1_02': -4.5741,
    'J2-1_03': -4.4376,
    'J2-1_04': -4.2209,
    'J2-1_05': -4.0976,
    'J2-1_06': -3.8808,
    'J2-1_07': -3.1619,
    'J2-1_08': -2.9453,
    'J2-1_09': -2.3469,
    'J2-1_10': -1.9290,
    'J2-1_11': -1.5888,
    'J2-1_12': -1.5516,
    'J2-1_13': -1.4523,
    'J2-1_14': -1.1465,
    'J2-1_15': -0.8065,
    'J2-1_16': -0.6532,
    'J2-1_17': -0.4694,
    'J2-1_18': -0.1767,
    'J2-1_19': 0.0000,
    'J2-1_20': 0.0071,
    'J2-1_21': 0.1137,
    'J2-1_22': 0.1291,
    'J2-1_23': 0.1617,
    'J2-1_24': 0.2239,
    'J2-1_25': 0.5237,
    'J2-1_26': 0.6384,
    'J2-1_27': 0.7405,
    'J2-1_28': 2.1394,
    'J2-1_29': 2.5158,
    'J2-1_30': 2.5444,
    'J2-1_31': 2.6225,
    'J2-1_32': 2.8844,
    'J2-1_33': 3.0325,
    'J2-1_34': 3.0990,
    'J2-1_35': 3.2981,
    'J2-1_36': 3.5091,
    'J2-1_37': 3.8148,
    'J2-1_38': 3.8201,
    'J2-1_39': 6.9891,
    'J2-1_40': 7.5057,
    ####### J 3-2
    'J3-2_01': -3.0666,
    'J3-2_02': -2.9296,
    'J3-2_03': -2.7221,
    'J3-2_04': -2.6563,
    'J3-2_05': -2.5270,
    'J3-2_06': -2.4010,
    'J3-2_07': -2.2535,
    'J3-2_08': -2.1825,
    'J3-2_09': -2.1277,
    'J3-2_10': -1.5862,
    'J3-2_11': -1.0158,
    'J3-2_12': -0.6131,
    'J3-2_13': -0.6093,
    'J3-2_14': -0.5902,
    'J3-2_15': -0.4872,
    'J3-2_16': -0.4725,
    'J3-2_17': -0.2757,
    'J3-2_18': -0.0697,
    'J3-2_19': -0.0616,
    'J3-2_20': -0.0022,
    'J3-2_21': 0.0000,
    'J3-2_22': 0.0143,
    'J3-2_23': 0.0542,
    'J3-2_24': 0.0561,
    'J3-2_25': 0.0575,
    'J3-2_26': 0.0687,
    'J3-2_27': 0.1887,
    'J3-2_28': 0.2411,
    'J3-2_29': 0.3781,
    'J3-2_30': 0.4620,
    'J3-2_31': 0.4798,
    'J3-2_32': 0.5110,
    'J3-2_33': 0.5540,
    'J3-2_34': 0.7808,
    'J3-2_35': 0.9066,
    'J3-2_36': 1.6382,
    'J3-2_37': 1.6980,
    'J3-2_38': 2.1025,
    'J3-2_39': 2.1236,
    'J3-2_40': 2.1815,
    'J3-2_41': 2.5281,
    'J3-2_42': 2.6458,
    'J3-2_43': 2.8052,
    'J3-2_44': 3.0320,
    'J3-2_45': 3.4963,
    }

line_strength_dict = {
    ####### J 1-0
    'J1-0_01': 0.025957,
    'J1-0_02': 0.065372,
    'J1-0_03': 0.019779,
    'J1-0_04': 0.004376,
    'J1-0_05': 0.034890,
    'J1-0_06': 0.071844,
    'J1-0_07': 0.259259,
    'J1-0_08': 0.156480,
    'J1-0_09': 0.028705,
    'J1-0_10': 0.041361,
    'J1-0_11': 0.013309,
    'J1-0_12': 0.056442,
    'J1-0_13': 0.156482,
    'J1-0_14': 0.028705,
    'J1-0_15': 0.037038,
    ####### J 2-1
    'J2-1_01': 0.008272,
    'J2-1_02': 0.005898,
    'J2-1_03': 0.031247,
    'J2-1_04': 0.013863,
    'J2-1_05': 0.013357,
    'J2-1_06': 0.010419,
    'J2-1_07': 0.000218,
    'J2-1_08': 0.000682,
    'J2-1_09': 0.000152,
    'J2-1_10': 0.001229,
    'J2-1_11': 0.000950,
    'J2-1_12': 0.000875,
    'J2-1_13': 0.002527,
    'J2-1_14': 0.000365,
    'J2-1_15': 0.000164,
    'J2-1_16': 0.021264,
    'J2-1_17': 0.031139,
    'J2-1_18': 0.000576,
    'J2-1_19': 0.200000,
    'J2-1_20': 0.001013,
    'J2-1_21': 0.111589,
    'J2-1_22': 0.088126,
    'J2-1_23': 0.142604,
    'J2-1_24': 0.011520,
    'J2-1_25': 0.027608,
    'J2-1_26': 0.012800,
    'J2-1_27': 0.066354,
    'J2-1_28': 0.013075,
    'J2-1_29': 0.003198,
    'J2-1_30': 0.061880,
    'J2-1_31': 0.004914,
    'J2-1_32': 0.035879,
    'J2-1_33': 0.011026,
    'J2-1_34': 0.039052,
    'J2-1_35': 0.019767,
    'J2-1_36': 0.004305,
    'J2-1_37': 0.001814,
    'J2-1_38': 0.000245,
    'J2-1_39': 0.000029,
    'J2-1_40': 0.000004,
    ####### J 3-2
    'J3-2_01': 0.001845,
    'J3-2_02': 0.001818,
    'J3-2_03': 0.003539,
    'J3-2_04': 0.014062,
    'J3-2_05': 0.011432,
    'J3-2_06': 0.000089,
    'J3-2_07': 0.002204,
    'J3-2_08': 0.002161,
    'J3-2_09': 0.000061,
    'J3-2_10': 0.000059,
    'J3-2_11': 0.000212,
    'J3-2_12': 0.000255,
    'J3-2_13': 0.000247,
    'J3-2_14': 0.000436,
    'J3-2_15': 0.010208,
    'J3-2_16': 0.000073,
    'J3-2_17': 0.007447,
    'J3-2_18': 0.000000,
    'J3-2_19': 0.000155,
    'J3-2_20': 0.000274,
    'J3-2_21': 0.174603,
    'J3-2_22': 0.018683,
    'J3-2_23': 0.135607,
    'J3-2_24': 0.100527,
    'J3-2_25': 0.124866,
    'J3-2_26': 0.060966,
    'J3-2_27': 0.088480,
    'J3-2_28': 0.001083,
    'J3-2_29': 0.094510,
    'J3-2_30': 0.014029,
    'J3-2_31': 0.007191,
    'J3-2_32': 0.022222,
    'J3-2_33': 0.047915,
    'J3-2_34': 0.015398,
    'J3-2_35': 0.000071,
    'J3-2_36': 0.000794,
    'J3-2_37': 0.001372,
    'J3-2_38': 0.007107,
    'J3-2_39': 0.016618,
    'J3-2_40': 0.009776,
    'J3-2_41': 0.000997,
    'J3-2_42': 0.000487,
    'J3-2_43': 0.000069,
    'J3-2_44': 0.000039,
    'J3-2_45': 0.000010,
}

# Get frequency dictionary in Hz based on the offset velocity and rest frequency
conv_J10=u.doppler_radio(freq_dict_cen['J1-0']*u.Hz)
conv_J21=u.doppler_radio(freq_dict_cen['J2-1']*u.Hz)
conv_J32=u.doppler_radio(freq_dict_cen['J3-2']*u.Hz)
freq_dict = {
    name: ((voff_lines_dict[name]*u.km/u.s).to(u.Hz, equivalencies=conv_J10).value) for name in voff_lines_dict.keys() if "J1-0" in name
    }
freq_dict.update({
    name: ((voff_lines_dict[name]*u.km/u.s).to(u.Hz, equivalencies=conv_J21).value) for name in voff_lines_dict.keys() if "J2-1" in name
    })
freq_dict.update({
    name: ((voff_lines_dict[name]*u.km/u.s).to(u.Hz, equivalencies=conv_J32).value) for name in voff_lines_dict.keys() if "J3-2" in name
    })

# relative_strength_total_degeneracy is not used in the CLASS implementation
# of the hfs fit. It is the sum of the degeneracy values for all hyperfines
# for a given line; it gives the relative weights between lines.
# Hyperfine weights are treated as normalized within one rotational transition.
w10 = sum(val for name,val in line_strength_dict.items() if 'J1-0' in name)
w21 = sum(val for name,val in line_strength_dict.items() if 'J2-1' in name)
w32 = sum(val for name,val in line_strength_dict.items() if 'J3-2' in name)
relative_strength_total_degeneracy = {
    name : w10 for name  in line_strength_dict.keys() if "J1-0" in name
    }
relative_strength_total_degeneracy.update({
    name : w21 for name  in line_strength_dict.keys() if "J2-1" in name
    })
relative_strength_total_degeneracy.update({
    name : w32 for name  in line_strength_dict.keys() if "J3-2" in name
    })

# Get the list of line names from the previous lists
line_names = [name for name in voff_lines_dict.keys()]

n2hp_vtau = hyperfine.hyperfinemodel(line_names, voff_lines_dict, freq_dict,
                                     line_strength_dict,
                                     relative_strength_total_degeneracy)
n2hp_vtau_fitter = n2hp_vtau.fitter
n2hp_vtau_vheight_fitter = n2hp_vtau.vheight_fitter
n2hp_vtau_tbg_fitter = n2hp_vtau.background_fitter

# RADEX part from old file

def n2hp_radex(xarr,
               density=4,
               column=13,
               xoff_v=0.0,
               width=1.0,
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
        print("density %20.12g column %20.12g: tau %20.12g tex %20.12g" % (density, column, tau, tex))

    if debug:
        import pdb; pdb.set_trace()

    return n2hp_vtau(xarr,Tex=tex,tau=tau,xoff_v=xoff_v,width=width,**kwargs)
