"""
===========
DCO+ fitter
===========
Reference for line params:

Dore (priv. comm.) line frequencies and strengths,
available at CDMS and Splatalogue

Rest freq Caselli, Dore 2005

"""
from . import hyperfine
import astropy.units as u

# line_names = ['J1-0', 'J2-1', 'J3-2',]
# line_names = ['J2-1', 'J3-2',]

freq_dict_cen ={
                'J1-0':  72039.3031e6,
                'J2-1': 144077.2804e6,
                'J3-2': 216112.5767e6,
               }

voff_lines_dict={
    ####### J 1-0
    'J1-0_01': -0.1977,
    'J1-0_02': 0.0000,
    'J1-0_03': 0.2568,
    ####### J 2-1
    'J2-1_01': -0.1086,
    'J2-1_02': -0.0901,
    'J2-1_03': -0.0098,
    'J2-1_04': 0.0000,
    'J2-1_05': 0.0385,
    'J2-1_06': 0.1373,
    ####### J 3-2
    'J3-2_01': -0.0771,
    'J3-2_02': -0.0171,
    'J3-2_03': -0.0046,
    'J3-2_04': 0.0000,
    'J3-2_05': 0.0086,
    'J3-2_06': 0.0810,
    }

line_strength_dict = {
    ####### J 1-0
    'J1-0_01': 0.3333,
    'J1-0_02': 0.5556,
    'J1-0_03': 0.1111,
    ####### J 2-1
    'J2-1_01': 0.083333,
    'J2-1_02': 0.111111,
    'J2-1_03': 0.250000,
    'J2-1_04': 0.466667,
    'J2-1_05': 0.005556,
    'J2-1_06': 0.083333,
    ####### J 3-2
    'J3-2_01': 0.0370,
    'J3-2_02': 0.2000,
    'J3-2_03': 0.2963,
    'J3-2_04': 0.4286,
    'J3-2_05': 0.0011,
    'J3-2_06': 0.0370,
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
relative_strength_total_degeneracy.update({
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

dcop_vtau = hyperfine.hyperfinemodel(line_names, voff_lines_dict, freq_dict,
                                     line_strength_dict,
                                     relative_strength_total_degeneracy)
dcop_vtau_fitter = dcop_vtau.fitter
dcop_vtau_vheight_fitter = dcop_vtau.vheight_fitter
dcop_vtau_tbg_fitter = dcop_vtau.background_fitter
