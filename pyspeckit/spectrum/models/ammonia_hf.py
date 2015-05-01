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
from . import model
from ammonia_constants import (line_names, freq_dict, aval_dict, ortho_dict,
                               voff_lines_dict, tau_wts_dict, line_labels)

from astropy import constants
from astropy import units as u
ckms = constants.c.to(u.km/u.s).value

relative_strength_total_degeneracy = collections.defaultdict(lambda: 1)

# sanity check:
for linename in line_names:
    assert len(voff_lines_dict[linename]) == len(tau_wts_dict[linename])

nh3_vtau = {linename:
            hyperfine.hyperfinemodel({lineid:lineid for lineid,name in
                                      enumerate(voff_lines_dict[linename])},
                                     {lineid:voff for lineid,voff in
                                      enumerate(voff_lines_dict[linename])},
                                     {lineid:freq_dict[linename]*(1-voff/ckms)
                                      for lineid,voff in
                                      enumerate(voff_lines_dict[linename])},
                                     {lineid:tauwt for lineid,tauwt in
                                      enumerate(tau_wts_dict[linename])},
                                     {lineid:1 for lineid,voff in
                                      enumerate(voff_lines_dict[linename])},
                                    )
            for linename in line_names}

def nh3_vtau_multimodel_generator(linenames):
    nlines = len(linenames)

    def nh3_vtau_multimodel(xarr, velocity, width, *args):
        assert len(args) == nlines*2
        models = [nh3_vtau[linename].hyperfine(xarr, Tex=tex, tau=tau,
                                               xoff_v=velocity, width=width)
                  for linename,tex,tau in zip(linenames, args[::2], args[1::2])]
        return np.array(models).sum(axis=0)


    mod = model.SpectralModel(nh3_vtau_multimodel, 2+nlines*2,
            parnames=['center','width'] + [x
                                           for ln in linenames
                                           for x in ('tex{0}'.format(ln),
                                                     'tau{0}'.format(ln))
                                          ],
            parlimited=[(False,False), (True,False),] + [(True, False),]*2*nlines, 
            parlimits=[(0,0), ]*(2+2*nlines),
            shortvarnames=["v","\\sigma",] + [x 
                                              for ln in linenames
                                              for x in
                                              ('T_{{ex}}({0})'.format(line_labels[ln]),
                                               '\\tau({0})'.format(line_labels[ln]))
                                             ],
            fitunits='Hz')

    return mod

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
