import pkgutil
import os
import importlib
from ...units import SpectroscopicAxis
import numpy as np
from ... import models
import pytest

models_path = os.path.dirname(models.__file__)
names = [name for _, name, _ in pkgutil.iter_modules([models_path])]

# this test doesn't work... it imports too much
@pytest.mark.xfail
@pytest.mark.parametrize('name',names)
def test_moments(name):
    """
    Try importing all the models and comparing their moments to their inputs
    """
    xarr = np.linspace(-100,100,200)
    rawdata = np.random.randn(xarr.size)

    model_name = 'pyspeckit.spectrum.models.' + name
    model_module = importlib.import_module(model_name)
    if hasattr(model_module, name+"_model"):
        model = getattr(model_module, name+'_model')
        model_instance = model()
        params = model_instance.moments(xarr, rawdata)
        print('params:', params)
        assert len(params) == model_instance.npars
        xarr = SpectroscopicAxis(xarr, unit='GHz')
        # if moments() returns wrong number of parameters this should fail
        print(model_instance.n_modelfunc(pars=params)(xarr))
