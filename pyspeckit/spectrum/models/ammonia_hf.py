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
import collections

from . import hyperfine
from . import radex_modelgrid
from ammonia_constants import (line_names, freq_dict, aval_dict, ortho_dict,
                               voff_lines_dict, tau_wts_dict)

relative_strength_total_degeneracy = collections.defaultdict(lambda: 1)

nh3_vtau = {linename:
            hyperfine.hyperfinemodel({lineid:name for lineid,name in
                                      enumerate(line_names[linename])},
                                     {lineid:voff for lineid,voff in
                                      enumerate(voff_lines_dict[linename])},
                                     {lineid:freq for lineid,freq in
                                      enumerate(freq_dict[linename])},
                                     {lineid:tauwt for lineid,tauwt in
                                      enumerate(tau_wts_dict[linename])},
                                     {lineid:1 for lineid,freq in
                                      enumerate(freq_dict[linename])},
                                    )
            for linename in line_names}

def nh3_vtau_multimodel_generator(linenames):
    def nh3_vtau_multimodel(xarr, velocity, width, *args):
        assert len(args) == len(linenames)*2
        models = [nh3_vtau[linename](xarr, tex, tau, velocity, width)

# nh3_vtau_fitter = nh3_vtau.fitter
# nh3_vtau_vheight_fitter = nh3_vtau.vheight_fitter

#nh3_radex = radex_modelgrid.radex_model(xarr, density=4, column=13, xoff_v=0.0,
#                                        width=1.0, grid_vwidth=1.0,
#                                        grid_vwidth_scale=False, texgrid=None,
#                                        taugrid=None, hdr=None,
#                                        path_to_texgrid='', path_to_taugrid='',
#                                        temperature_gridnumber=3,
#                                        modelfunc=nh3_vtau, debug=False,
#                                        verbose=False)
