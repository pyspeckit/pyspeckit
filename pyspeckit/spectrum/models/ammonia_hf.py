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
from pyspeckit.mpfit import mpfit
import fitter
import matplotlib.cbook as mpcb
import copy

from ammonia_constants import (line_names, freq_dict, aval_dict, ortho_dict,
                               voff_lines_dict, tau_wts_dict)


nh3_vtau = hyperfine.hyperfinemodel(line_names, voff_lines_dict, freq_dict,
                                    line_strength_dict,
                                    relative_strength_total_degeneracy)
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
