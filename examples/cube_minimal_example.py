import pyspeckit
from pyspeckit.cubes.tests.test_cubetools import make_test_cube
import matplotlib.pylab as plt
import numpy as np

# generate a test spectral cube (10x10, with a 100 spectral channels)
make_test_cube((100,10,10), outfile='test.fits')
spc = pyspeckit.Cube('test.fits')

# do a crude noise estimate on the 30 edge channels
rmsmap = np.vstack([spc.cube[:15], spc.cube[85:]]).std(axis=0)

# get a cube of moments
spc.momenteach(vheight=False)

# fit each pixel taking its moment as an initial guess
spc.fiteach(fittype = 'gaussian',
            guesses = spc.momentcube,
            errmap = rmsmap,
            signal_cut = 3, # ignore pixels with SNR<3
            blank_value = np.nan,
            start_from_point=(5,5))
spc.mapplot()
# show the fitted amplitude
spc.show_fit_param(0, cmap='viridis')
plt.show()
