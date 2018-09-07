from __future__ import print_function
import pyspeckit
import numpy as np
from astropy import units as u

from pyspeckit.spectrum.models import ammonia

xarr = np.linspace(-40, 40, 300) * u.km/u.s
oneonemod = ammonia.ammonia(xarr.to(u.GHz, u.doppler_radio(ammonia.freq_dict['oneone']*u.Hz)),)
twotwomod = ammonia.ammonia(xarr.to(u.GHz, u.doppler_radio(ammonia.freq_dict['twotwo']*u.Hz)),)

sp11 = pyspeckit.Spectrum(xarr=xarr, data=oneonemod, unit=u.K,
                          xarrkwargs={'refX': ammonia.freq_dict['oneone']*u.Hz},
                          header={})
sp22 = pyspeckit.Spectrum(xarr=xarr, data=twotwomod, unit=u.K,
                          xarrkwargs={'refX': ammonia.freq_dict['twotwo']*u.Hz},
                          header={})


input_dict={'oneone':sp11, 'twotwo':sp22,}
spf, specout = pyspeckit.wrappers.fitnh3.fitnh3tkin(input_dict, dobaseline=False)
print(specout.specfit.modelpars)
print(specout.specfit.parinfo)

spf2, specout2 = pyspeckit.wrappers.fitnh3.fitnh3tkin(input_dict,
                                                      dobaseline=True,
                                                      baselinekwargs={'exclude':[-30,30]*u.km/u.s})
print(specout.specfit.modelpars)
print(specout.specfit.parinfo)
