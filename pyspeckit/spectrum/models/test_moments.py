import pkgutil
import os
import importlib
from ..models.inherited_gaussfitter import gaussian_fitter
from ..units import SpectroscopicAxis
import numpy as np
from . import models

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
		print model_instance.n_modelfunc(pars=params)(xarr)
	except Exception as e:
		print e
