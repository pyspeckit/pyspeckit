from astropy import constants
from astropy import units as u
import numpy as np

TCMB = 2.7326 # K

num_to_name = {0: 'zero',
               1: 'one',
               2: 'two',
               3: 'three',
               4: 'four',
               5: 'five',
               6: 'six',
               7: 'seven',
               8: 'eight',
               9: 'nine'}

line_names = ['oneone','twotwo','threethree','fourfour', "fivefive", "sixsix",
              "sevenseven", "eighteight", 'ninenine']

# indices of the level population array for ortho/para lines
line_name_indices = {'oneone': 0,
                     'twotwo': 1,
                     'fourfour': 2,
                     'fivefive': 3,
                     'sevenseven': 4,
                     'eighteight': 5,
                     'tenten': 6,
                     'eleveneleven': 7,
                     'zerozero': 0,
                     'threethree': 1,
                     'sixsix': 2,
                     'ninenine': 3,
                    }

              #'ninenine']
line_labels = {'oneone': '1-1',
               'twotwo': '2-2',
               'threethree': '3-3',
               'fourfour': '4-4',
               'fivefive': '5-5',
               'sixsix': '6-6',
               'sevenseven': '7-7',
               'eighteight': '8-8',
               'ninenine': '9-9',
              }

freq_dict = { 
    'oneone':     23.6944955e9, # Issue 91: Erik's custom frequency
    'twotwo':     23.722633335e9,
    'threethree': 23.8701296e9,
    'fourfour':   24.1394169e9,
    # fivefive - eighteight come from TopModel on splatalogue
    # see ammonia_offset_calculation.py
    'fivefive': 24.53299e9,
    'sixsix': 25.05603e9,
    'sevenseven': 25.71518e9,
    'eighteight': 26.51898e9,
    'ninenine': 27.477943e9,
    }
aval_dict = {
    'oneone':     1.712e-7,  #64*!pi**4/(3*h*c**3)*nu11**3*mu0**2*(1/2.)
    'twotwo':     2.291e-7,  #64*!pi**4/(3*h*c**3)*nu22**3*mu0**2*(2/3.)
    'threethree': 2.625e-7,  #64*!pi**4/(3*h*c**3)*nu33**3*mu0**2*(3/4.)
    'fourfour':   3.167e-7,  #64*!pi**4/(3*h*c**3)*nu44**3*mu0**2*(4/5.)
    'fivefive': 3.099109e-07,
    'sixsix': 3.395797e-07,
    'sevenseven': 3.747560e-07,
    'eighteight': 4.175308e-07,
    'ninenine': 2.602045278086488e-07,
    }
ortho_dict = {
    'oneone':     False,
    'twotwo':     False,
    'threethree': True,
    'fourfour':   False,
    'fivefive':   False,
    'sixsix':   True,
    'sevenseven':   False,
    'eighteight':   False,
    'ninenine':   True,
    }
#n_ortho = np.arange(0,28,3) # 0..3..27
#n_para = np.array([x for x in range(28) if x % 3 != 0])

voff_lines_dict = {
    'oneone': [19.8513, 19.3159, 7.88669, 7.46967, 7.35132, 0.460409, 0.322042,
        -0.0751680, -0.213003, 0.311034, 0.192266, -0.132382, -0.250923, -7.23349,
        -7.37280, -7.81526, -19.4117, -19.5500],
    'twotwo':[26.5263, 26.0111, 25.9505, 16.3917, 16.3793, 15.8642, 0.562503,
        0.528408, 0.523745, 0.0132820, -0.00379100, -0.0132820, -0.501831,
        -0.531340, -0.589080, -15.8547, -16.3698, -16.3822, -25.9505, -26.0111,
        -26.5263],
    'threethree':[29.195098, 29.044147, 28.941877, 28.911408, 21.234827,
        21.214619, 21.136387, 21.087456, 1.005122, 0.806082, 0.778062,
        0.628569, 0.016754, -0.005589, -0.013401, -0.639734, -0.744554,
        -1.031924, -21.125222, -21.203441, -21.223649, -21.076291, -28.908067,
        -28.938523, -29.040794, -29.191744],
    'fourfour':[  0.        , -30.49783692,  30.49783692,   0., 24.25907811,
        -24.25907811,   0.        ],
    # fivefive - eighteight come from TopModel on splatalogue
    # see ammonia_offset_calculation.py
    'fivefive': [31.4053287863, 26.0285409785, 0.0, 0.0, 0.0, -25.9063412556, -31.2831290633],
    'sixsix': [31.5872901302, 27.0406347326, 0.0, 0.0, 0.0, -26.9209859064, -31.4676413039],
    'sevenseven': [31.3605314845, 27.3967468359, 0.0, 0.0, 0.0, -27.5133287373, -31.477113386],
    'eighteight': [30.9752235915, 27.4707274918, 0.0, 0.0, 0.0, -27.5837757531, -30.9752235915],
    'ninenine': [0],
}

tau_wts_dict = {
    'oneone': [0.0740740, 0.148148, 0.0925930, 0.166667, 0.0185190, 0.0370370,
        0.0185190, 0.0185190, 0.0925930, 0.0333330, 0.300000, 0.466667,
        0.0333330, 0.0925930, 0.0185190, 0.166667, 0.0740740, 0.148148],
    'twotwo': [0.00418600, 0.0376740, 0.0209300, 0.0372090, 0.0260470,
        0.00186000, 0.0209300, 0.0116280, 0.0106310, 0.267442, 0.499668,
        0.146512, 0.0116280, 0.0106310, 0.0209300, 0.00186000, 0.0260470,
        0.0372090, 0.0209300, 0.0376740, 0.00418600],
    'threethree': [0.012263, 0.008409, 0.003434, 0.005494, 0.006652, 0.008852,
        0.004967, 0.011589, 0.019228, 0.010387, 0.010820, 0.009482, 0.293302,
        0.459109, 0.177372, 0.009482, 0.010820, 0.019228, 0.004967, 0.008852,
        0.006652, 0.011589, 0.005494, 0.003434, 0.008409, 0.012263],
    # fivefive - eighteight come from TopModel on splatalogue
    # see ammonia_offset_calculation.py
    'fourfour': [0.2431, 0.0162, 0.0162, 0.3008, 0.0163, 0.0163, 0.3911],
    'fivefive': [0.0109080940831, 0.0109433143618, 0.311493418617, 0.261847767275, 0.382955997218, 0.0109433143618, 0.0109080940831],
    'sixsix': [0.0078350431801, 0.00784948916416, 0.317644539734, 0.274246689798, 0.376739705779, 0.00784948916416, 0.0078350431801],
    'sevenseven': [0.00589524944656, 0.00590204051181, 0.371879455317, 0.321515700951, 0.283010263815, 0.00590204051181, 0.00589524944656],
    'eighteight': [0.00459516014524, 0.00459939439378, 0.324116135075, 0.289534720829, 0.367960035019, 0.00459939439378, 0.00459516014524],
    'ninenine': [1],
}

#ckms = 2.99792458e5
ckms = constants.c.to(u.km/u.s).value
ccms = constants.c.to(u.cm/u.s).value
#Degeneracies
# g1 = 1
# g2 = 1
#h = 6.6260693e-27
h = constants.h.cgs.value
#kb = 1.3806505e-16
kb = constants.k_B.cgs.value
# Dipole Moment in cgs (1.476 Debeye)
#mu0 = 1.476e-18

# Generate Partition Functions
nlevs = 51
jv=np.arange(nlevs)
ortho = jv % 3 == 0
para = ~ortho
Jpara = jv[para]
Jortho = jv[ortho]
Brot = 298117.06e6
Crot = 186726.36e6
