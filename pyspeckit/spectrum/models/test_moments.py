# this script 1: loads all models in its directory,
# and runs a loop where for each model 2: it imports it,
# 3: instantiates it (using the convention that the model has a class named NAME-OF-MODEL_model)
# then 4: calls its moments() method with some fake data
# and finally 5: calls the models .n_modelfunc() method
# using the params from the moments() call.
# If the call does not raise an Exception, then the moments() method
# returns the correct number of parameters.

import pkgutil
import os
import importlib
from pyspeckit.spectrum.units import SpectroscopicAxis
import numpy as np
import pyspeckit.spectrum.models as models
import pytest

models_path = os.path.dirname(models.__file__)
xarr = np.linspace(-100,100,200)
names = [name for _, name, _ in pkgutil.iter_modules([models_path])]
for name in names:
    try:
        model_name = 'pyspeckit.spectrum.models.' + name
        model_module = importlib.import_module(model_name)
        model = getattr(model_module, name+'_model')
        model_instance = model()
        moments = getattr(model, 'moments')
        rawdata = np.random.randn(xarr.size)
        params = model_instance.moments(xarr, rawdata)
        print 'params:', params
        xarr = SpectroscopicAxis(xarr)
        # if moments() returns wrong number of parameters this should fail 
        model_instance.n_modelfunc(pars=params)(xarr)
    except Exception as e:
        pytest.fail(e)
