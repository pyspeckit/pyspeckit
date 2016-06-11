"""
===========
NH2D fitter: ortho- and para- in the same file, but not modeled together
===========
Reference for line params:

Shah & Wooten (2001) line frequencies and line strengths.

TODO: Obtaine up-to-date values from CDMS and Splatalogue

Ronak Y. Shah and Alwyn Wootten ApJ 554, 933-947 (2001)
http://adsabs.harvard.edu/abs/2001ApJ...554..933S

"""
from . import hyperfine
import astropy.units as u

freq_dict_cen ={
                'o-1_01-1_11':  85.926263e9,
                'p-1_01-1_11': 110.153599e9,
                'o-1_01-0_00': 332.82251e9,
                'p-1_01-0_00': 332.78189e9,
               }

freq_dict={
    ####### ortho-NH2D J=1_01-1_11
    'o-1_01-1_11_01': 85.9247829e9,
    'o-1_01-1_11_02': 85.9257031e9,
    'o-1_01-1_11_03': 85.9262703e9,
    'o-1_01-1_11_04': 85.9263165e9,
    'o-1_01-1_11_05': 85.9268837e9,
    'o-1_01-1_11_06': 85.9277345e9,
    ####### ortho-NH2D J=1_01-0_00
    'o-1_01-1_11_01': 332.7809447e9,
    'o-1_01-1_11_02': 332.7817955e9,
    'o-1_01-1_11_03': 332.7823627e9,
    ####### para-NH2D J=1_01-1_11
    'p-1_01-0_00_01': 110.152084e9,
    'p-1_01-0_00_02': 110.152995e9,
    'p-1_01-0_00_03': 110.153599e9,
    'p-1_01-0_00_04': 110.153599e9,
    'p-1_01-0_00_05': 110.154222e9,
    'p-1_01-0_00_06': 110.155053e9,
    ####### para-NH2D J=1_01-0_00
    'p-1_01-0_00_01': 332.8215595e9,
    'p-1_01-0_00_02': 332.8224149e9,
    'p-1_01-0_00_03': 332.8229853e9,
    }

line_strength_dict = {
    ####### ortho-NH2D J=1_01-1_11
    'o-1_01-1_11_01': 0.111,
    'o-1_01-1_11_02': 0.139,
    'o-1_01-1_11_03': 0.417,
    'o-1_01-1_11_04': 0.083,
    'o-1_01-1_11_05': 0.139,
    'o-1_01-1_11_06': 0.111,
    ####### ortho-NH2D J=1_01-0_00
    'o-1_01-0_00_01': 0.111,
    'o-1_01-0_00_02': 0.556,
    'o-1_01-0_00_03': 0.333,
    ####### para-NH2D J=1_01-1_11
    'p-1_01-1_11_01': 0.111,
    'p-1_01-1_11_02': 0.139,
    'p-1_01-1_11_03': 0.417,
    'p-1_01-1_11_04': 0.083,
    'p-1_01-1_11_05': 0.139,
    'p-1_01-1_11_06': 0.111,
    ####### para-NH2D J=1_01-0_00
    'p-1_01-0_00_01': 0.111,
    'p-1_01-0_00_02': 0.556,
    'p-1_01-0_00_03': 0.333,
}

# freq_dict = {
#     'J2-1': (voff_lines_dict['J2-1']*u.km/u.s).to(u.GHz, equivalencies=u.doppler_radio(freq_dict_cen['J2-1']*u.Hz)).value,
#     'J3-2': (voff_lines_dict['J3-2']*u.km/u.s).to(u.GHz, equivalencies=u.doppler_radio(freq_dict_cen['J3-2']*u.Hz)).value,
# }
# Get offset velocity dictionary in km/s based on the lines frequencies and rest frequency
conv_o1_1=u.doppler_radio(freq_dict_cen['o-1_01-1_11']*u.Hz)
conv_p1_1=u.doppler_radio(freq_dict_cen['p-1_01-1_11']*u.Hz)
conv_o1_0=u.doppler_radio(freq_dict_cen['o-1_01-0_00']*u.Hz)
conv_p1_0=u.doppler_radio(freq_dict_cen['p-1_01-0_00']*u.Hz)

voff_lines_dict = {
    name: ((freq_dict[name]*u.Hz).to(u.km/u.s, equivalencies=conv_o1_1).value) for name in freq_dict.keys() if "o-1_01-1_11" in name
    }
voff_lines_dict.update({
	name: ((freq_dict[name]*u.Hz).to(u.km/u.s, equivalencies=conv_p1_1).value) for name in freq_dict.keys() if "p-1_01-1_11" in name
    })
voff_lines_dict.update({
	name: ((freq_dict[name]*u.Hz).to(u.km/u.s, equivalencies=conv_p1_1).value) for name in freq_dict.keys() if "o-1_01-0_00" in name
    })
voff_lines_dict.update({
	name: ((freq_dict[name]*u.Hz).to(u.km/u.s, equivalencies=conv_p1_1).value) for name in freq_dict.keys() if "p-1_01-0_00" in name
    })

# I don't know yet how to use this parameter... in CLASS it does not exist
# Note to Jaime: this is the sum of the degeneracy values for all hyperfines
# for a given line; it gives the relative weights between the J=2-1 and J=3-2
# lines, for example (the hyperfine weights are treated as normalized within
# one rotational transition)
wo1_1 = sum(val for name,val in voff_lines_dict.items() if 'o-1_01-1_11' in name)
wp1_1 = sum(val for name,val in voff_lines_dict.items() if 'p-1_01-1_11' in name)
wo1_0 = sum(val for name,val in voff_lines_dict.items() if 'o-1_01-0_00' in name)
wp1_0 = sum(val for name,val in voff_lines_dict.items() if 'p-1_01-0_00' in name)
relative_strength_total_degeneracy = {
    name : wo1_1 for name  in voff_lines_dict.keys() if "o-1_01-1_11" in name
    }
relative_strength_total_degeneracy.update({
    name : wp1_1 for name  in voff_lines_dict.keys() if "p-1_01-1_11" in name
    })
relative_strength_total_degeneracy.update({
    name : wo1_0 for name  in voff_lines_dict.keys() if "o-1_01-0_00" in name
    })
relative_strength_total_degeneracy.update({
    name : wp1_0 for name  in voff_lines_dict.keys() if "p-1_01-0_00" in name
    })
# Get the list of line names from the previous lists
line_names = [name for name in voff_lines_dict.keys()]


nh2d_vtau = hyperfine.hyperfinemodel(line_names, voff_lines_dict, freq_dict,
                                     line_strength_dict,
                                     relative_strength_total_degeneracy)
nh2d_vtau_fitter = nh2d_vtau.fitter
nh2d_vtau_vheight_fitter = nh2d_vtau.vheight_fitter
nh2d_vtau_tbg_fitter = nh2d_vtau.background_fitter
