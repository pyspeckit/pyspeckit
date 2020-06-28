import pyspeckit as psk
from pyspeckit.spectrum import models
from astropy.table import Table
from spectral_cube import SpectralCube
import numpy as np
import matplotlib.pyplot as plt
import despotic
import pyspeckit.spectrum.readers.read_class
import os
import shutil


if not os.path.exists('ph2cogrid.fits'):
    if not os.path.exists('protostellarCore.desp'):
        despotic_install_path = (os.path.split(despotic.__file__))[0]
        shutil.copy(despotic_install_path+'/cloudfiles/protostellarCore.desp',os.getcwd())
        models.formaldehyde_mm.build_despotic_grids(gridfile='ph2cogrid.fits', DvUpper=10)

t = Table.read('ph2cogrid.fits')

# This returns interpolating functions that take physical parameters
# and returns values for Tex, Tau for the three mm transitions.
f1, f2, f3 = models.formaldehyde_mm.formaldehyde_mm_despotic_functions(t)

# Instantiate that fitter!
formaldehyde_fitter=models.model.SpectralModel(models.formaldehyde_mm.formaldehyde_mm_despotic,
                                               5, parnames=['temperature', 'column', 'density',
                                                            'center', 'width'],
                                               parvalues=[50,12,5.0,0,2],
                                               parlimited=[(True, True), (True, True),
                                                           (True, True), (False, False),
                                                           (True, False)],
                                               parlimits=[(5,205), (10,17),
                                                          (2,7), (0,0), (0,0)],
                                               parsteps=[0.01, 0.01, 0.1, 0, 0],
                                               fitunits='Hz',
                                               h2co_303_202=f1, # interpolation of (Tex, tau)
                                               h2co_322_221=f2,
                                               h2co_321_220=f3,
                                               shortvarnames=("T", "N", "n", "v", "\\sigma"))

sp = pyspeckit.readers.read_class.class_to_spectra('example_h2co_mm_spectrum.apex')
sp.data *= 1/0.75 # T_A* -> T_MB
sp.unit = "$T_{MB}$"
# estimate the error from the data
# sp.error[:] = sp.stats((2.183e2,2.184e2))['std']

sp.Registry.add_fitter('formaldehyde_mm_despotic', formaldehyde_fitter, 5)
#plot fit for all 3 ('both')
sp.plotter(figure=1)
sp.specfit(fittype='formaldehyde_mm_despotic',
           guesses=[95, 14.5, 4, 0.0, 4.0],
           limits=[(10,300), (11,15), (2,7), (-20,150), (1, 10)],
           limited=[(True, True)]*5,
           fixed=[False, False, True, False, False])

sp.plotter.savefig('test_fitting_figure_01.png')
