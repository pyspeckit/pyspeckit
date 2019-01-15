import pyspeckit
import os
from pyspeckit.spectrum.models import nh2d
import numpy as np

import astropy.units as u

if not os.path.exists('o-nh2d_spec.fits'):
    import astropy.utils.data as aud
    from astropy.io import fits
    f = aud.download_file('https://github.com/pyspeckit/pyspeckit-example-files/raw/master/o-nh2d_spec.fits')
    with fits.open(f) as ff:
        ff.writeto('o-nh2d_spec.fits')

# Load the spectrum 
spec = pyspeckit.Spectrum('o-nh2d_spec.fits')
# Determine rms from line free section and load into cube
rms = np.std(spec.data[10:370])
spec.error[:] = rms
# setup spectral axis
spec.xarr.refX = 85.92627*u.GHz
spec.xarr.velocity_convention = 'radio'
spec.xarr.convert_to_unit('km/s')
# define useful shortcuts for True and False
F=False
T=True
# Setup of matplotlib
import matplotlib.pyplot as plt
plt.ion()
# Add NH2D fitter
spec.Registry.add_fitter('nh2d_vtau', pyspeckit.models.nh2d.nh2d_vtau_fitter,4)
# run spectral fit using some reasonable guesses
spec.specfit(fittype='nh2d_vtau', guesses=[6.51, 4.4, 0.214, 0.1088], 
    verbose_level=4, signal_cut=1.5, limitedmax=[F,T,T,T], limitedmin=[T,T,T,T], 
    minpars=[0, 0, -1, 0.05], maxpars=[30.,50.,1,0.5], fixed=[F,F,F,F])
# plot best fit
spec.plotter(errstyle='fill')
spec.specfit.plot_fit()
#save figure
plt.savefig('example_o-NH2D.png')
