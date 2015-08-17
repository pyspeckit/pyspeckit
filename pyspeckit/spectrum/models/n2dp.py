"""
===========
N2D+ fitter
===========
Reference for line params: 

L. Dore, P. Caselli, S. Beninati, T. Bourke, P. C. Myers and G. Cazzoli A&A 413, 1177-1181 (2004) 

http://adsabs.harvard.edu/abs/2004A%26A...413.1177D

"""
import numpy as np
from pyspeckit.mpfit import mpfit
from .. import units
from . import fitter,model,modelgrid
import matplotlib.cbook as mpcb
import copy
import hyperfine
from pyspeckit.specwarnings import warn
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

# line_names = ['oneone_f10', 'oneone_f01', 'oneone_f22', 'oneone_f21',
#               'oneone_f12', 'oneone_f11', 'twotwo_f11', 'twotwo_f12',
#               'twotwo_f21', 'twotwo_f32', 'twotwo_f33', 'twotwo_f22',
#               'twotwo_f23']

# Name is JN1-N2 - F'1,F1 -> F'2,F2
freq_dict ={
'J2-1_22_21':154214.845e6, 
'J2-1_22_23':154215.040e6, 
'J2-1_21_21':154215.168e6, 
'J2-1_23_23':154215.288e6, 
'J2-1_22_22':154215.331e6, 
'J2-1_11_01':154215.425e6, 
'J2-1_23_22':154215.579e6, 
'J2-1_12_01':154215.628e6, 
'J2-1_21_22':154215.653e6, 
'J2-1_10_01':154215.888e6,
'J2-1_22_11':154216.754e6, 
'J2-1_33_23':154216.818e6, 
'J2-1_22_12':154216.889e6, 
'J2-1_21_11':154217.076e6, 
'J2-1_33_22':154217.109e6, 
'J2-1_32_21':154217.131e6, 
'J2-1_23_12':154217.137e6, 
'J2-1_34_23':154217.206e6, 
'J2-1_21_10':154217.502e6, 
'J2-1_32_22':154217.617e6, 
'J2-1_12_23':154218.120e6, 
'J2-1_11_11':154219.631e6, 
'J2-1_11_12':154219.766e6, 
'J2-1_12_11':154219.834e6, 
'J2-1_12_12':154219.969e6, 
'J2-1_11_10':154220.058e6, 
'J2-1_10_11':154220.094e6, 
'J3-2_32_32':231319.945e6, 
'J3-2_34_34':231319.995e6, 
'J3-2_33_33':231320.017e6, 
'J3-2_22_12':231321.253e6, 
'J3-2_22_11':231321.456e6, 
'J3-2_21_10':231321.499e6, 
'J3-2_44_34':231321.530e6, 
'J3-2_33_23':231321.547e6, 
'J3-2_23_12':231321.617e6, 
'J3-2_33_22':231321.795e6, 
'J3-2_32_21':231321.908e6, 
'J3-2_34_23':231321.914e6, 
'J3-2_44_33':231321.918e6, 
'J3-2_43_32':231321.919e6, 
'J3-2_21_11':231321.961e6,
'J3-2_45_34':231321.966e6, 
'J3-2_32_22':231322.231e6, 
'J3-2_43_33':231322.426e6, 
'J3-2_22_21':231324.012e6, 
'J3-2_22_23':231324.086e6,
'J3-2_22_22':231324.334e6,
'J3-2_23_23':231324.450e6,
'J3-2_21_21':231324.517e6, 
'J3-2_23_22':231324.698e6,
'J3-2_21_22':231324.839e6}

line_strength_dict = { # effectively the degeneracy per rotation state...
'J2-1_22_21':0.002, 
'J2-1_22_23':0.004, 
'J2-1_21_21':0.020, 
'J2-1_23_23':0.039, 
'J2-1_22_22':0.011, 
'J2-1_11_01':0.036, 
'J2-1_23_22':0.005, 
'J2-1_12_01':0.062, 
'J2-1_21_22':0.003, 
'J2-1_10_01':0.013,
'J2-1_22_11':0.066, 
'J2-1_33_23':0.013, 
'J2-1_22_12':0.027, 
'J2-1_21_11':0.012, 
'J2-1_33_22':0.143, 
'J2-1_32_21':0.088, 
'J2-1_23_12':0.112, 
'J2-1_34_23':0.200, 
'J2-1_21_10':0.031, 
'J2-1_32_22':0.021, 
'J2-1_12_23':0.003, 
'J2-1_11_11':0.010, 
'J2-1_11_12':0.013, 
'J2-1_12_11':0.014, 
'J2-1_12_12':0.031, 
'J2-1_11_10':0.006, 
'J2-1_10_11':0.008, 
'J3-2_32_32':0.010, 
'J3-2_34_34':0.017, 
'J3-2_33_33':0.007, 
'J3-2_22_12':0.015, 
'J3-2_22_11':0.048, 
'J3-2_21_10':0.022, 
'J3-2_44_34':0.007, 
'J3-2_33_23':0.014, 
'J3-2_23_12':0.094, 
'J3-2_33_22':0.089, 
'J3-2_32_21':0.061, 
'J3-2_34_23':0.125, 
'J3-2_44_33':0.136, 
'J3-2_43_32':0.101, 
'J3-2_21_11':0.019,
'J3-2_45_34':0.175, 
'J3-2_32_22':0.007, 
'J3-2_43_33':0.010, 
'J3-2_22_21':0.002, 
'J3-2_22_23':0.002,
'J3-2_22_22':0.011,
'J3-2_23_23':0.014,
'J3-2_21_21':0.004, 
'J3-2_23_22':0.002,
'J3-2_21_22':0.002}

line_names = freq_dict.keys()

# http://adsabs.harvard.edu/abs/1971ApJ...169..429T has the most accurate freqs
# http://adsabs.harvard.edu/abs/1972ApJ...174..463T [twotwo]
central_freq_dict = { 
    'J2-1': 154217.1805e6,
    'J3-2': 231321.9119e6
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
        'twotwo_f33': 41.48,
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
        'twotwo_f11':14.48847881e9 - 19.97e3,
        'twotwo_f12':14.48847881e9 -  7.03e3,
        'twotwo_f21':14.48847881e9 -  2.20e3,
        'twotwo_f32':14.48847881e9 +  0.12e3,
        'twotwo_f33':14.48847881e9 +  0.89e3,
        'twotwo_f22':14.48847881e9 + 10.74e3,
        'twotwo_f23':14.48847881e9 + 11.51e3,
        'threethree_f22':28.97478e9,
        'threethree_f44':28.97480e9,
        'threethree_f33':28.97481e9,
        }
freq_dict = dict(hf_freq_dict.items() + central_freq_dict.items())
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
