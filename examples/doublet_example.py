import pyspeckit

# Read in J000002.09+155254.1 spectrum, a nice emission-line galaxy
sp = pyspeckit.Spectrum('../tests/SIIdoublet.fits')

SIIa = sp.speclines.optical.lines['SIIa'][0]
SIIb = sp.speclines.optical.lines['SIIb'][0]

offset = SIIb - SIIa

# Let's have a look.
sp.plotter()

raw_input('\nLet\'s do a simple continuum subtraction (continue)')

sp.baseline(subtract = False)

raw_input('\nLet\'s zoom in on the SII doublet (continue)')

sp.baseline()
sp.plotter.savefig('doublet_example_fullspectrum.png')

z = 0.02        # I already know the redshift, so lets take that into account

# Zoom in on SII doublet
sp.plotter(xmin = SIIa * (1 + z) - 75, xmax = SIIb * (1 + z) + 75, ymin = -10, ymax = 60)

# Guess positions to be rest wavelengths times factor (1 + z)
guesses = [100, SIIa * (1 + z), 5, 100, SIIb * (1 + z), 5]
tied = ['', '', '', '', 'p[1] + %g' % offset, '']

sp.specfit(guesses = guesses, tied = tied)
sp.plotter.savefig('doublet_example_SII.png')
              
raw_input('\nHooray! The doublet has been fit. ')

