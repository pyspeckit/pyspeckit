import spectrum

if not 'interactive' in globals():
    interactive=False

sp = spectrum.Spectrum('sample_13CO.fits')

print "Does it have an axis? ",sp.plotter.axis
sp.plotter()
print "How about now? ",sp.plotter.axis


# set the baseline to zero to prevent variable-height fitting
# (if you don't do this, the best fit to the spectrum is dominated by the
# background level)
sp.baseline.order = 0
sp.specfit()
sp.plotter.figure.savefig('fits_gaussfit.png')
print "Guesses: ", sp.specfit.guesses
print "Best fit: ", sp.specfit.modelpars

# # Fit a baseline, excluding the velocities with data, and don't subtract it
# sp.baseline(exclude=[12000,98000],subtract=False)
# print "Baseline: ",sp.baseline.baselinepars
# print "Excludepix: ",sp.baseline.excludepix
# print "EQW: ",sp.specfit.EQW()

sp.specfit(interactive=True)
if interactive: raw_input('Press enter to print guesses and best fit and end code')
sp.plotter.figure.savefig('fits_interactive_fit.png')
print "Guesses: ", sp.specfit.guesses
print "Best fit: ", sp.specfit.modelpars

# don't try this for a zero-baseline spectrum print "EQW: ",sp.specfit.EQW()

print "Attempting to write to test.fits: "
sp.write('test.fits')

sptest = spectrum.Spectrum('test.fits')
sptest.plotter()
if interactive: raw_input('Test fits done.')

print "Attempting to write to test.hdf5: "
sp.write('test.hdf5')

#from matplotlib import pyplot

