"""
===========
N2D+ fitter
===========
Reference for line params:

Dore (priv. comm.) line frequencies in CDMS,
line strength can also be obtained from Splatalogue

L. Dore, P. Caselli, S. Beninati, T. Bourke, P. C. Myers and G. Cazzoli A&A 413, 1177-1181 (2004)
http://adsabs.harvard.edu/abs/2004A%26A...413.1177D

L. Pagani, F. Daniel, and M. L. Dubernet A\%A 494, 719-727 (2009)
DOI: 10.1051/0004-6361:200810570

In this version added N2D+(1-0). Data taken from the articles above.

Module API
^^^^^^^^^^
"""
from . import hyperfine
import astropy.units as u

# line_names = ['J1-0', 'J2-1', 'J3-2',]
# line_names = ['J2-1', 'J3-2',]

freq_dict_cen ={
                'J1-0':  77109.6162e6,
                'J2-1': 154217.1805e6,
                'J3-2': 231321.9119e6,
               }

voff_lines_dict={
    ####### J 1-0
    'J1-0_01': -9.6730,
    'J1-0_02': -9.6730,
    'J1-0_03': -9.6730,
    'J1-0_04': -0.7537,
    'J1-0_05': -0.7537,
    'J1-0_06': -0.7537,
    'J1-0_07': 0.0000,
    'J1-0_08': 1.1314,
    'J1-0_09': 1.1314,
    'J1-0_10': 6.6519,
    'J1-0_11': 6.6519,
    'J1-0_12': 6.6519,
    'J1-0_13': 7.1766,
    'J1-0_14': 7.1766,
    'J1-0_15': 8.3080,
    ####### J 2-1
    'J2-1_01': -5.6031,
    'J2-1_02': -5.5332,
    'J2-1_03': -5.3617,
    'J2-1_04': -5.0993,
    'J2-1_05': -4.9677,
    'J2-1_06': -4.7052,
    'J2-1_07': -3.8195,
    'J2-1_08': -3.5571,
    'J2-1_09': -2.8342,
    'J2-1_10': -2.3388,
    'J2-1_11': -1.9449,
    'J2-1_12': -1.9002,
    'J2-1_13': -1.7733,
    'J2-1_14': -1.3965,
    'J2-1_15': -1.0025,
    'J2-1_16': -0.7968,
    'J2-1_17': -0.5740,
    'J2-1_18': -0.2311,
    'J2-1_19': -0.0085,
    'J2-1_20': 0.0000,
    'J2-1_21': 0.1351,
    'J2-1_22': 0.1457,
    'J2-1_23': 0.1886,
    'J2-1_24': 0.2538,
    'J2-1_25': 0.6165,
    'J2-1_26': 0.7541,
    'J2-1_27': 0.8789,
    'J2-1_28': 2.5594,
    'J2-1_29': 3.0143,
    'J2-1_30': 3.0632,
    'J2-1_31': 3.1579,
    'J2-1_32': 3.4572,
    'J2-1_33': 3.6394,
    'J2-1_34': 3.7234,
    'J2-1_35': 3.9567,
    'J2-1_36': 4.2049,
    'J2-1_37': 4.5817,
    'J2-1_38': 4.6054,
    'J2-1_39': 8.4164,
    'J2-1_40': 9.0414,
    ####### J 3-2
    'J3-2_01': -3.7164,
    'J3-2_02': -3.5339,
    'J3-2_03': -3.2997,
    'J3-2_04': -3.2130,
    'J3-2_05': -3.0633,
    'J3-2_06': -2.8958,
    'J3-2_07': -2.7424,
    'J3-2_08': -2.6466,
    'J3-2_09': -2.5748,
    'J3-2_10': -1.9177,
    'J3-2_11': -1.2333,
    'J3-2_02': -0.7628,
    'J3-2_13': -0.7590,
    'J3-2_14': -0.7306,
    'J3-2_15': -0.5953,
    'J3-2_16': -0.5765,
    'J3-2_17': -0.3419,
    'J3-2_18': -0.0925,
    'J3-2_19': -0.0210,
    'J3-2_20': 0.0000,
    'J3-2_21': 0.0065,
    'J3-2_22': 0.0616,
    'J3-2_23': 0.0618,
    'J3-2_24': 0.0675,
    'J3-2_25': 0.0748,
    'J3-2_26': 0.2212,
    'J3-2_27': 0.2691,
    'J3-2_28': 0.4515,
    'J3-2_29': 0.5422,
    'J3-2_30': 0.5647,
    'J3-2_31': 0.6050,
    'J3-2_32': 0.6596,
    'J3-2_33': 0.9222,
    'J3-2_34': 1.0897,
    'J3-2_35': 1.9586,
    'J3-2_36': 2.0471,
    'J3-2_37': 2.5218,
    'J3-2_38': 2.5500,
    'J3-2_39': 2.6156,
    'J3-2_40': 3.0245,
    'J3-2_41': 3.1786,
    'J3-2_42': 3.3810,
    'J3-2_43': 3.6436,
    'J3-2_44': 4.2066,
    }

line_strength_dict = {
    ####### J 1-0
    'J1-0_01': 0.026018,
    'J1-0_02': 0.065408,
    'J1-0_03': 0.019683,
    'J1-0_04': 0.004390,
    'J1-0_05': 0.035006,
    'J1-0_06': 0.071714,
    'J1-0_07': 0.259259,
    'J1-0_08': 0.156212,
    'J1-0_09': 0.028973,
    'J1-0_10': 0.041311,
    'J1-0_11': 0.013379,
    'J1-0_12': 0.056422,
    'J1-0_13': 0.156214,
    'J1-0_14': 0.028973,
    'J1-0_15': 0.037038,
    ####### J 2-1
    'J2-1_01': 0.008262,
    'J2-1_02': 0.005907,
    'J2-1_03': 0.031334,
    'J2-1_04': 0.013833,
    'J2-1_05': 0.013341,
    'J2-1_06': 0.010384,
    'J2-1_07': 0.000213,
    'J2-1_08': 0.000675,
    'J2-1_09': 0.000150,
    'J2-1_10': 0.001202,
    'J2-1_11': 0.000963,
    'J2-1_12': 0.000878,
    'J2-1_13': 0.002533,
    'J2-1_14': 0.000362,
    'J2-1_15': 0.000162,
    'J2-1_16': 0.021268,
    'J2-1_17': 0.031130,
    'J2-1_18': 0.000578,
    'J2-1_19': 0.001008,
    'J2-1_20': 0.200000,
    'J2-1_21': 0.111666,
    'J2-1_22': 0.088138,
    'J2-1_23': 0.142511,
    'J2-1_24': 0.011550,
    'J2-1_25': 0.027472,
    'J2-1_26': 0.012894,
    'J2-1_27': 0.066406,
    'J2-1_28': 0.013082,
    'J2-1_29': 0.003207,
    'J2-1_30': 0.061847,
    'J2-1_31': 0.004932,
    'J2-1_32': 0.035910,
    'J2-1_33': 0.011102,
    'J2-1_34': 0.038958,
    'J2-1_35': 0.019743,
    'J2-1_36': 0.004297,
    'J2-1_37': 0.001830,
    'J2-1_38': 0.000240,
    'J2-1_39': 0.000029,
    'J2-1_40': 0.000004,
    ####### J 3-2
    'J3-2_01': 0.001842,
    'J3-2_02': 0.001819,
    'J3-2_03': 0.003544,
    'J3-2_04': 0.014100,
    'J3-2_05': 0.011404,
    'J3-2_06': 0.000088,
    'J3-2_07': 0.002201,
    'J3-2_08': 0.002153,
    'J3-2_09': 0.000059,
    'J3-2_10': 0.000058,
    'J3-2_11': 0.000203,
    'J3-2_12': 0.000259,
    'J3-2_13': 0.000248,
    'J3-2_14': 0.000437,
    'J3-2_15': 0.010215,
    'J3-2_16': 0.000073,
    'J3-2_17': 0.007445,
    'J3-2_18': 0.000155,
    'J3-2_19': 0.000272,
    'J3-2_20': 0.174603,
    'J3-2_21': 0.018678,
    'J3-2_22': 0.100524,
    'J3-2_23': 0.135563,
    'J3-2_24': 0.124910,
    'J3-2_25': 0.060970,
    'J3-2_26': 0.088513,
    'J3-2_27': 0.001085,
    'J3-2_28': 0.094480,
    'J3-2_29': 0.013955,
    'J3-2_30': 0.007236,
    'J3-2_31': 0.022222,
    'J3-2_32': 0.047921,
    'J3-2_33': 0.015427,
    'J3-2_34': 0.000070,
    'J3-2_35': 0.000796,
    'J3-2_36': 0.001373,
    'J3-2_37': 0.007147,
    'J3-2_38': 0.016574,
    'J3-2_39': 0.009776,
    'J3-2_40': 0.000995,
    'J3-2_41': 0.000491,
    'J3-2_42': 0.000067,
    'J3-2_43': 0.000039,
    'J3-2_44': 0.000010,
}

# freq_dict = {
#     'J2-1': (voff_lines_dict['J2-1']*u.km/u.s).to(u.GHz, equivalencies=u.doppler_radio(freq_dict_cen['J2-1']*u.Hz)).value,
#     'J3-2': (voff_lines_dict['J3-2']*u.km/u.s).to(u.GHz, equivalencies=u.doppler_radio(freq_dict_cen['J3-2']*u.Hz)).value,
# }

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

# I don't know yet how to use this parameter... in CLASS it does not exist
# Note to Jaime: this is the sum of the degeneracy values for all hyperfines
# for a given line; it gives the relative weights between the J=2-1 and J=3-2
# lines, for example (the hyperfine weights are treated as normalized within
# one rotational transition)
w10 = sum(val for name,val in line_strength_dict.items() if 'J1-0' in name)
w21 = sum(val for name,val in line_strength_dict.items() if 'J2-1' in name)
w32 = sum(val for name,val in line_strength_dict.items() if 'J3-2' in name)
relative_strength_total_degeneracy = {
    name : w10 for name  in line_strength_dict.keys() if "J1-0" in name
    }
relative_strength_total_degeneracy.update ({
    name : w21 for name  in line_strength_dict.keys() if "J2-1" in name
    })
relative_strength_total_degeneracy.update({
    name : w32 for name  in line_strength_dict.keys() if "J3-2" in name
    })

# Get the list of line names from the previous lists
line_names = [name for name in voff_lines_dict.keys()]

#     'J2-1': np.array([1]*len(voff_lines_dict['J2-1'])),
#     'J3-2': np.array([1]*len(voff_lines_dict['J3-2'])),
# }
# aval_dict = {
#     # 'J1-0': 10**(-4.90770),
#     'J2-1': 10**(-3.92220),
#     'J3-2': 10**(-3.35866),
# }

n2dp_vtau = hyperfine.hyperfinemodel(line_names, voff_lines_dict, freq_dict,
                                     line_strength_dict,
                                     relative_strength_total_degeneracy)
n2dp_vtau_fitter = n2dp_vtau.fitter
n2dp_vtau_vheight_fitter = n2dp_vtau.vheight_fitter
n2dp_vtau_tbg_fitter = n2dp_vtau.background_fitter
