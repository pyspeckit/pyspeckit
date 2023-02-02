from astropy import constants
from astropy import units as u
import numpy as np

# References for spectroscopy:
# cdms catalogue https://cdms.astro.uni-koeln.de/cgi-bin/cdmsinfo?file=e018503.cat
# In particular, for frequencies: J. T. Hougen, 1972, J. Chem. Phys. 57, 4207;
# The rotational constants are average values. They were derived from the data set described in:
# P. Chen, J. C. Pearson, H. M. Pickett, S. Matsuura, and G.A. Blake, 2006, J. Mol. Spectrosc. 236, 116.

# The full hyperfine structure is based on calculations performed by Luca Bizzocchi <luca.bizzocchi@unibo.it>, performed from a complete reanalysis
# of the literature molecular beam data (Kukolich 1967, 1968; Hougen 1972)

num_to_name = {0: 'zero',
               1: 'one',
               2: 'two'}


line_names = ['oneone', 'twotwo']

# indices of the level population array for ortho/para lines
line_name_indices = {'oneone': 0,
                     'twotwo': 1,
                    }


line_labels = {'oneone': '1-1',
               'twotwo': '2-2'}

freq_dict = {
    'oneone': 22624.9295e6,
    'twotwo': 22649.8434e6
}

aval_dict = {
    'oneone': 10**-6.8355,
    'twotwo': 10**-6.7093
}

ortho_dict = {
    'oneone': False,
    'twotwo': False
}

voff_lines_dict = {
    'oneone': [0.556880, 0.535679, 0.432325, 0.397874, 0.294520, 0.139488,
               -0.028793, -0.048669, -0.153348, -0.230202, -0.250077,
               -0.354757, -0.444860],
    'twotwo': [0.6578, 0.6115, 0.3759, 0.2806, 0.2356, 0.0000, -0.2356,
               -0.2806, -0.3759, -0.6115, -0.6578]
}


tau_wts_dict = {
    'oneone': [0.086485, 0.093960, 0.000386, 0.040499, 0.042829, 0.069190,
               0.017691, 0.064914, 0.153229, 0.312494, 0.050622, 0.053548,
               0.014153],
    'twotwo': [0.029846, 0.008717, 0.001239, 0.009040, 0.028779, 0.844579,
               0.028779, 0.009040, 0.001239, 0.008717, 0.029846]
}


ckms = constants.c.to(u.km/u.s).value
ccms = constants.c.to(u.cm/u.s).value
h = constants.h.cgs.value
kb = constants.k_B.cgs.value


# Generate Partition Functions
nlevs = 51  ##TO DO: is this correct??
jv = np.arange(nlevs)
ortho = jv % 3 == 0
para = ~ortho
Jpara = jv[para]
Jortho = jv[ortho]

Brot = 297390.97e6
Crot = 186734.0e6
