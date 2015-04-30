import numpy as np
import pyspeckit
from astropy import units as u
from pyspeckit.spectrum.models import ammonia_constants, ammonia, ammonia_hf

# Generate a synthetic spectrum based off of 3 NH3 lines
xarr11 = pyspeckit.units.SpectroscopicAxis(np.linspace(-30, 30, 100)*u.km/u.s,
                                           velocity_convention='radio',
                                           refX=ammonia_constants.freq_dict['oneone']).as_unit(u.GHz)
xarr22 = pyspeckit.units.SpectroscopicAxis(np.linspace(-30, 30, 100)*u.km/u.s,
                                           velocity_convention='radio',
                                           refX=ammonia_constants.freq_dict['twotwo']).as_unit(u.GHz)
xarr33 = pyspeckit.units.SpectroscopicAxis(np.linspace(-30, 30, 100)*u.km/u.s,
                                           velocity_convention='radio',
                                           refX=ammonia_constants.freq_dict['threethree']).as_unit(u.GHz)
xarr = pyspeckit.units.SpectroscopicAxes([xarr11,xarr22,xarr33])

synthspec = (ammonia.ammonia(xarr, tkin=20, ntot=15, fortho=0.5) +
             ammonia.ammonia(xarr, tkin=50, ntot=14, fortho=0.5))

spectrum = pyspeckit.Spectrum(xarr=xarr, data=synthspec)

spectrum.plotter()

spectrum.specfit.Registry.add_fitter('nh3_vtau', ammonia_hf.nh3_vtau.fitter, ammonia_hf.nh3_vtau.fitter.npars)
# parnames=['Tex','tau','center','width'],
spectrum.specfit(fittype='nh3_vtau', guesses=[5, 0.5, 0.5, 0.5])
