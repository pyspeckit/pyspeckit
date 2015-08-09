"""
simple cube fitting example
"""
from __future__ import print_function

import numpy as np
import pyspeckit
from spectral_cube import SpectralCube
from astropy import wcs
from astropy import units as u

# Create a new WCS object so we can instantiate the SpectralCube
mywcs = wcs.WCS(naxis=3)

# Set up a tangent projection
# You would normally read this from a file!!
mywcs.wcs.crpix = [1,2,21.0]
mywcs.wcs.cdelt = np.array([-0.066667, 0.066667, 500])
mywcs.wcs.crval = [290.9250, 14.5092, 60000]
mywcs.wcs.ctype = ["RA---TAN", "DEC--TAN", 'VELO']
mywcs.wcs.cunit = ['deg', 'deg', 'm/s']

# Create a synthetic X-dimension in km/s
xarr = np.linspace(50, 70, 41) # km/s

# Define a line width, which will vary across our image
# It will increase from 1 km/s to 4 km/s over the X-direction (RA)
sigma = np.outer(np.linspace(1,1.5,2), np.ones(4)).T

# Define a line center, which will vary in the opposite direction,
# along increasing Y-direction (declination)
centroid = np.outer(np.ones(2), np.linspace(58, 62, 4)).T

data = np.exp(-(np.tile(xarr, (2, 4, 1)).T - centroid)**2 / (2.*sigma**2))
cube = SpectralCube(data=data, wcs=mywcs)

# Sanity checks: do the moments accurately recover the inputs?
assert (np.abs(cube.moment1().to(u.km/u.s).value - centroid).max()) < 1e-5
assert (np.abs(cube.moment2().to(u.km**2/u.s**2).value - sigma**2).max()) < 1e-5

# Create a pyspeckit cube
pcube = pyspeckit.Cube(cube=cube)

# For convenience, convert the X-axis to km/s
# (WCSLIB automatically converts to m/s even if you give it km/s)
pcube.xarr.convert_to_unit(u.km/u.s)

# Set up the fitter by doing a preliminary fit
pcube.specfit(fittype='gaussian', guesses='moments')

# Fit each spectrum with a gaussian
# First, assemble the guesses:
guesses = np.array([cube.max(axis=0).value,
                    cube.moment1(axis=0).to(u.km/u.s).value,
                    (cube.moment2(axis=0)**0.5).to(u.km/u.s).value])
# (the second moment is in m^2/s^2, but we want km/s

# Do the fit!
pcube.fiteach(guesses=guesses, # pass in the guess array
              # tell it where to start the fitting (center pixel in this case)
              start_from_point=(5,5),
              # Paralellize the fits?
              multicore=4,
              fittype='gaussian',
             )

# Then you can access the fits via parcube:
assert np.all(pcube.parcube[0,:,:] == 1)
assert np.all(pcube.parcube[1,:,:] == centroid)
assert np.all(pcube.parcube[2,:,:] == sigma)

# Can also fit non-parallelized
pcube.fiteach(guesses=guesses, # pass in the guess array
              # tell it where to start the fitting (center pixel in this case)
              start_from_point=(5,5),
              # Paralellize the fits?
              multicore=1,
              fittype='gaussian',
             )

assert np.all(pcube.parcube[0,:,:] == 1)
assert np.all(pcube.parcube[1,:,:] == centroid)
assert np.all(pcube.parcube[2,:,:] == sigma)
