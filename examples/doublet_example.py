from __future__ import print_function
import pyspeckit

# Read in J000002.09+155254.1 spectrum, a nice emission-line galaxy
sp = pyspeckit.Spectrum('SIIdoublet.fits')

# Read in rest wavelengths of SII lines.  If you didn't know the names already, 
# you could do sp.speclines.optical.lines.keys() to see what is available.
SIIa = sp.speclines.optical.lines['SIIa'][0]
SIIb = sp.speclines.optical.lines['SIIb'][0]

# Wavelength difference between doublet lines - use to tie positions together
offset = SIIb - SIIa

# Let's have a look at the spectrum
sp.plotter()

# Let's do a simple continuum subtraction (continue)

# Plot the baseline fit
sp.baseline(subtract = False)

# Let's zoom in on the SII doublet (continue)

# Subtract the baseline fit and save
sp.baseline(subtract = True)
sp.plotter.savefig('doublet_example_fullspectrum.png')

# Guess for the redshift - otherwise we'll end up with the Halpha-NII complex
z = 0.02        

# Zoom in on SII doublet
sp.plotter(xmin = SIIa * (1 + z) - 75, xmax = SIIb * (1 + z) + 75, ymin = -10, ymax = 60)

# Guess amplitudes to be 100, positions to be rest wavelengths 
# times factor (1 + z), and widths to be 5 Angstroms
guesses = [100, SIIa * (1 + z), 5, 100, SIIb * (1 + z), 5]
tied = ['', '', '', '', 'p[1] + %g' % offset, '']

# Do the fit, and plot it
sp.specfit(guesses = guesses, tied = tied, quiet = False)
sp.plotter.savefig('doublet_example_SII.png')
              
# Hooray! The doublet has been fit.

SIIb_obs = sp.specfit.modelpars[-2]

print('Our guess for the redshift was z = 0.02.')
print('The redshift, as derived by the line shift, is z = %g' % ((SIIb_obs / SIIb) - 1))
