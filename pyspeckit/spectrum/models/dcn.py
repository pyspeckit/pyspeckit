"""
===========
DCN fitter:
===========
Reference for line params:

Line strength and frequencies taken from CDMS in Aug. 17th 2023:
DCN v=0, tag: 28509

"""
from . import hyperfine
import astropy.units as u

freq_dict_cen = {
'''
Frequency from CDMS without HFS
'''
                'J1-0': 72414.6936e6,
                'J2-1': 144828.0015e6,
                'J3-2': 217238.5378e6,
                }

freq_dict = {
    ####### J 1-0
    'J1-0_01': 72413.5040e6,
    'J1-0_02': 72414.9330e6,
    'J1-0_03': 72417.0280e6,
    ####### J 2-1
    'J2-1_01': 144826.5757e6, 
    'J2-1_02': 144826.8216e6,
    'J2-1_03': 144828.0011e6,
    'J2-1_04': 144828.1106e6,
    'J2-1_05': 144828.9072e6,
    'J2-1_06': 144830.3326e6,
    ####### J 3-2
    'J3-2_01': 217236.9990e6, 
    'J3-2_02': 217238.3000e6,
    'J3-2_03': 217238.5550e6,
    'J3-2_04': 217238.6120e6,
    'J3-2_05': 217239.0790e6,
    'J3-2_06': 217240.6220e6,
    }

line_strength_dict = {
    ####### J 1-0
    'J1-0_01': 0.33331025,
    'J1-0_02': 0.55558213,
    'J1-0_03': 0.11110761,
    ####### J 2-1
    'J2-1_01': 0.08333412,
    'J2-1_02': 0.11110224,
    'J2-1_03': 0.24999011,
    'J2-1_04': 0.46668391,
    'J2-1_05': 0.0055555,
    'J2-1_06': 0.08333412,
    ####### J 3-2
    'J3-2_01': 0.0370362, 
    'J3-2_02': 0.19999833,
    'J3-2_03': 0.2962964,
    'J3-2_04': 0.42857454,
    'J3-2_05': 0.00105834,
    'J3-2_06': 0.0370362,
    }


# Get offset velocity dictionary in km/s based on the lines frequencies and rest frequency
conv_J10 = u.doppler_radio(freq_dict_cen['J1-0']*u.Hz)
conv_J21 = u.doppler_radio(freq_dict_cen['J2-1']*u.Hz)
conv_J32 = u.doppler_radio(freq_dict_cen['J3-2']*u.Hz)

voff_lines_dict = {
    name: ((freq_dict[name]*u.Hz).to(u.km/u.s, equivalencies=conv_J10).value) for name in freq_dict.keys() if "J1-0" in name
    }
voff_lines_dict.update({
	name: ((freq_dict[name]*u.Hz).to(u.km/u.s, equivalencies=conv_J21).value) for name in freq_dict.keys() if "J2-1" in name
    })
voff_lines_dict.update({
	name: ((freq_dict[name]*u.Hz).to(u.km/u.s, equivalencies=conv_J32).value) for name in freq_dict.keys() if "J3-2" in name
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

dcn_vtau = hyperfine.hyperfinemodel(line_names, voff_lines_dict, freq_dict,
                                     line_strength_dict,
                                     relative_strength_total_degeneracy)
dcn_vtau_fitter = dcn_vtau.fitter
dcn_vtau_vheight_fitter = dcn_vtau.vheight_fitter
dcn_vtau_tbg_fitter = dcn_vtau.background_fitter
