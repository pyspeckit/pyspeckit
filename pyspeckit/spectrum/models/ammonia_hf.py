"""
===================================================
Ammonia inversion transition: Hyperfine-only fitter
===================================================

.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>

Module API
^^^^^^^^^^

"""
import numpy as np
import matplotlib.cbook as mpcb
import copy
import collections

from ...mpfit import mpfit
from . import fitter
from . import hyperfine
from . import radex_modelgrid
from . import model
from .ammonia_constants import (line_names, freq_dict, aval_dict, ortho_dict,
                                TCMB,
                                voff_lines_dict, tau_wts_dict, line_labels)

from astropy import constants
from astropy import units as u
ckms = constants.c.to(u.km/u.s).value

# sanity check:
for linename in line_names:
    assert len(voff_lines_dict[linename]) == len(tau_wts_dict[linename])

# For each individual inversion line, create a Hyperfine model
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
                                     {lineid:sum(tau_wts_dict[linename]) for lineid,voff in
                                      enumerate(voff_lines_dict[linename])},
                                    )
            for linename in line_names}

def nh3_vtau_multimodel_generator(linenames):
    """
    If you want to use multiple hyperfines for the same spectrum, use this
    generator.
    It is useful if you want N independent tau/tex values but the same velocity
    and linewidth

    Parameters
    ----------
    linenames : list
        A list of line names from the set ('oneone', ..., 'eighteight')

    Returns
    -------
    model : `model.SpectralModel`
        A SpectralModel class build from N different metastable inversion
        hyperfine models
    """
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
            parlimits=[(0,0), (0,0), ] + [(TCMB, 0), (0, 0)] * nlines,
            shortvarnames=["v","\\sigma",] + [x
                                              for ln in linenames
                                              for x in
                                              ('T_{{ex}}({0})'.format(line_labels[ln]),
                                               '\\tau({0})'.format(line_labels[ln]))
                                             ],
            fitunit='Hz')

    return mod
