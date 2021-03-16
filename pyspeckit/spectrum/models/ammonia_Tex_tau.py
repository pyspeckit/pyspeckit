"""
===========
ammonia fitter
===========
Reference for line params:

Dore (priv. comm.) line frequencies in CDMS,
line strength can also be obtained from Splatalogue

Petrashkevich and Punanova by analogy with other models.

Module API
^^^^^^^^^^
"""
from . import hyperfine
import astropy.units as u

# line_names = ['J1-0', 'J2-1', 'J3-2',]
# line_names = ['J2-1', 'J3-2',]

freq_dict_cen ={
                'J1-1':  23.6944955e9
               }

voff_lines_dict={
    ####### J 1-1
    'J1-1_01':19.8513,                 #-19.3329,
    'J1-1_02':19.3159,                 #-7.8915,
    'J1-1_03':7.88669,                 # -7.4742,
    'J1-1_04':7.46967,                 # -7.3551,
    'J1-1_05':7.35132,                 #-0.46257,
    'J1-1_06':0.460409,                 #-0.32334,
    'J1-1_07':0.322042,                 #-0.30885,
    'J1-1_08':0.311034,                 #-0.189624,
    'J1-1_09':0.192266,                 #0.074,
    'J1-1_10':-0.0751680,                #0.133131,
    'J1-1_11': -0.213003,               #0.2133,
    'J1-1_12': -0.132382,                #0.25236,
    'J1-1_13': -0.250923,                #7.2393,
    'J1-1_14': -7.23349,                #7.3788,
    'J1-1_15':-7.37280,               #7.8207,
    'J1-1_16': -7.81526,               #19.4226,
    'J1-1_17': -19.4117, 
    'J1-1_18': -19.5500,                      #19.5621,
}


line_strength_dict = {
    ####### J 1-0
    'J1-1_01':0.0740740, 
    'J1-1_02':0.148148,
    'J1-1_03':0.0925930, 
    'J1-1_04':0.166667,
    'J1-1_05':0.0185190, 
    'J1-1_06':0.0370370, 
    'J1-1_07':0.0185190, 
    'J1-1_08':0.0333330, 
    'J1-1_09':0.300000,
    'J1-1_10':0.0185190, 
    'J1-1_11':0.0925930, 
    'J1-1_12':0.466667, 
    'J1-1_13':0.0333330, 
    'J1-1_14':0.0925930,
    'J1-1_15':0.0185190, 
    'J1-1_16':0.166667, 
    'J1-1_17':0.0740740, 
    'J1-1_18':0.148148,
}

# freq_dict = {
#     'J2-1': (voff_lines_dict['J2-1']*u.km/u.s).to(u.GHz, equivalencies=u.doppler_radio(freq_dict_cen['J2-1']*u.Hz)).value,
#     'J3-2': (voff_lines_dict['J3-2']*u.km/u.s).to(u.GHz, equivalencies=u.doppler_radio(freq_dict_cen['J3-2']*u.Hz)).value,
# }

# Get frequency dictionary in Hz based on the offset velocity and rest frequency
conv_J10=u.doppler_radio(freq_dict_cen['J1-1']*u.Hz)
#conv_J21=u.doppler_radio(freq_dict_cen['J2-1']*u.Hz)
#conv_J32=u.doppler_radio(freq_dict_cen['J3-2']*u.Hz)
freq_dict = {
    name: ((voff_lines_dict[name]*u.km/u.s).to(u.Hz, equivalencies=conv_J10).value) for name in voff_lines_dict.keys() if "J1-1" in name
    }
'''
freq_dict.update({
    name: ((voff_lines_dict[name]*u.km/u.s).to(u.Hz, equivalencies=conv_J21).value) for name in voff_lines_dict.keys() if "J2-1" in name
    })
freq_dict.update({
    name: ((voff_lines_dict[name]*u.km/u.s).to(u.Hz, equivalencies=conv_J32).value) for name in voff_lines_dict.keys() if "J3-2" in name
    })
'''
# I don't know yet how to use this parameter... in CLASS it does not exist
# Note to Jaime: this is the sum of the degeneracy values for all hyperfines
# for a given line; it gives the relative weights between the J=2-1 and J=3-2
# lines, for example (the hyperfine weights are treated as normalized within
# one rotational transition)
w10 = sum(val for name,val in line_strength_dict.items() if 'J1-1' in name)
#w21 = sum(val for name,val in line_strength_dict.items() if 'J2-1' in name)
#w32 = sum(val for name,val in line_strength_dict.items() if 'J3-2' in name)
relative_strength_total_degeneracy = {
    name : w10 for name  in line_strength_dict.keys() if "J1-1" in name
    }
'''
relative_strength_total_degeneracy.update ({
    name : w21 for name  in line_strength_dict.keys() if "J2-1" in name
    })
relative_strength_total_degeneracy.update({
    name : w32 for name  in line_strength_dict.keys() if "J3-2" in name
    })
'''
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

amm_vtau = hyperfine.hyperfinemodel(line_names, voff_lines_dict, freq_dict,
                                     line_strength_dict,
                                     relative_strength_total_degeneracy)
amm_vtau_fitter = amm_vtau.fitter
amm_vtau_vheight_fitter = amm_vtau.vheight_fitter
amm_vtau_tbg_fitter = amm_vtau.background_fitter
