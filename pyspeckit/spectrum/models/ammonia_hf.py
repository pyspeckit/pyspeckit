"""
========================================
Ammonia inversion transition TKIN fitter
========================================

Ammonia inversion transition TKIN fitter translated from Erik Rosolowsky's
http://svn.ok.ubc.ca/svn/signals/nh3fit/

.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>

Module API
^^^^^^^^^^

"""
import numpy as np
from mpfit import mpfit
import fitter
import matplotlib.cbook as mpcb
import copy


freq_dict = { 
    'oneone':     23.694506e9,
    'twotwo':     23.722633335e9,
    'threethree': 23.8701296e9,
    'fourfour':   24.1394169e9,
    }
aval_dict = {
    'oneone':     1.712e-7,  #64*!pi**4/(3*h*c**3)*nu11**3*mu0**2*(1/2.)
    'twotwo':     2.291e-7,  #64*!pi**4/(3*h*c**3)*nu22**3*mu0**2*(2/3.)
    'threethree': 2.625e-7,  #64*!pi**4/(3*h*c**3)*nu33**3*mu0**2*(3/4.)
    'fourfour':   3.167e-7,  #64*!pi**4/(3*h*c**3)*nu44**3*mu0**2*(4/5.)
    }
ortho_dict = {
    'oneone':     False,
    'twotwo':     False,
    'threethree': True,
    'fourfour':   False,
    }
n_ortho = np.arange(0,28,3) # 0..3..27
n_para = np.array([x for x in range(28) if x % 3 != 0])

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
        -24.25907811,   0.        ]
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
    'fourfour': [0.2431, 0.0162, 0.0162, 0.3008, 0.0163, 0.0163, 0.3911]}


line_names = freq_dict.keys()
# line_names = ['oneone','twotwo','threethree','fourfour']

ckms = units.speedoflight_ms / 1e3 #2.99792458e5
# doesn't work for nh3voff_lines_dict = dict([(k,(v-93.176261e9)/93.176261e9*ckms) for k,v in freq_dict.iteritems()])

nh3_vtau = hyperfine.hyperfinemodel(line_names, voff_lines_dict, freq_dict, line_strength_dict)
nh3_vtau_fitter = nh3_vtau.fitter
nh3_vtau_vheight_fitter = nh3_vtau.vheight_fitter

nh3_radex = radex_modelgrid(xarr, density=4, column=13, xoff_v=0.0, width=1.0, 
        grid_vwidth=1.0,
        grid_vwidth_scale=False,
        texgrid=None,
        taugrid=None,
        hdr=None,
        path_to_texgrid='',
        path_to_taugrid='',
        temperature_gridnumber=3,
        modelfunc=nh3_vtau,
        debug=False,
        verbose=False)
