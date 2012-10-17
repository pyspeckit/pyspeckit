"""

test_measurements.py

Author: Jordan Mirocha
Affiliation: University of Colorado at Boulder
Created on: Sat Oct  6 14:24:37 2012

Description: Simple test of SDSS data, and the measurements class.

"""

import pyspeckit

# Open file
sp = pyspeckit.Spectrum('./sample_sdss.fits', specnum = 1, errspecnum = 2)

sp.plotter()
print 'Optical spectrum of SDSS J095751.40+105834.6.'

sp.plotter(xmin = 5000, xmax = 5300)
print 'Zoom-in on the OIII-Hbeta complex.'

# Let's look at the OIII-H-beta complex
sp.specfit.selectregion(5000, 5300)

guesses = [100, 5050, 5, 40, 5150, 5, 120, 5200, 5]
sp.specfit(guesses = guesses, negamp = False)

print 'Fit 3 Gaussians.\nCall the measurements class for automatic line identification and such.'
sp.measure(z = sp.header.get('Z'), fluxnorm = 1e-17, restframe = True)

print 'Line   Wavelength (A)  FWHM (A)  Flux (erg/s/cm^2/Ang)  Luminosity (erg/s)'
for line in sp.measurements.lines:
    print line, sp.measurements.lines[line]['pos'], sp.measurements.lines[line]['fwhm'], \
        sp.measurements.lines[line]['flux'], sp.measurements.lines[line]['lum']
