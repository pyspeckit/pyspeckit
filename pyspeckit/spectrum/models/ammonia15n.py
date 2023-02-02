"""
========================================
15NH3 inversion transition TROT fitter
========================================

Ammonia inversion transition TROT fitter translated from Erik Rosolowsky's
https://github.com/low-sky/nh3fit

.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>
modified for 15n-bearing isotopologue by:: Elena Redaelli <eredaelli@mpe.mpg.de>

Module API
^^^^^^^^^^

#############TODO################
Implement the equivalent of the cold_ammonia_model for 14NH3 in 15NH3
to allow for Tkin to be fitted (instead of Trot)

"""
from __future__ import division

import numpy as np

from .ammonia15n_constants import (line_names, freq_dict, aval_dict,
                                   ortho_dict, voff_lines_dict, tau_wts_dict)

from .ammonia import ammonia, ammonia_model

TCMB = 2.7315 # K

def ammonia15n(xarr, line_names=line_names, freq_dict=freq_dict,
               aval_dict=aval_dict, ortho_dict=ortho_dict,
               voff_lines_dict=voff_lines_dict, tau_wts_dict=tau_wts_dict,
               **kwargs):
    """
    See docstring for ammonia model
    """
    return ammonia(xarr, line_names=line_names, freq_dict=freq_dict,
                   aval_dict=aval_dict, ortho_dict=ortho_dict,
                   voff_lines_dict=voff_lines_dict, tau_wts_dict=tau_wts_dict,
                   **kwargs)


class ammonia15n_model(ammonia_model):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.modelfunc = ammonia15n


class ammonia15n_model_vtau(ammonia15n_model):
    def __init__(self,
                 parnames=['trot', 'tex', 'tau', 'width', 'xoff_v', 'fortho'],
                 **kwargs):
        super().__init__(parnames=parnames, **kwargs)


class ammonia15n_model_background(ammonia15n_model):
    def __init__(self, **kwargs):
        super().__init__(npars=7, parnames=['trot', 'tex', 'ntot', 'width',
                                            'xoff_v', 'fortho',
                                            'background_tb'])


class ammonia15n_model_restricted_tex(ammonia15n_model):
    def __init__(self,
                 parnames=['trot', 'tex', 'ntot', 'width', 'xoff_v', 'fortho',
                           'delta'],
                 ammonia15n_func=ammonia15n,
                 **kwargs):
        super().__init__(npars=7, parnames=parnames, **kwargs)

        def ammonia15n_dtex(*args, **kwargs):
            """
            Strip out the 'delta' keyword
            """
            # for py2 compatibility, must strip out manually
            delta = kwargs.pop('delta') if 'delta' in kwargs else None
            np.testing.assert_allclose(kwargs['trot'] - kwargs['tex'],
                                       delta)
            return ammonia15n_func(*args, **kwargs)
        self.modelfunc = ammonia15n_dtex
