import spectrum

sp = spectrum.Spectrum('sample_13CO.fits')

print "Does it have an axis? ",sp.plotter.axis
sp.plotter()
print "How about now? ",sp.plotter.axis


# set the baseline to zero to prevent variable-height fitting
# (if you don't do this, the best fit to the spectrum is dominated by the
# background level)
sp.baseline.order = 0
sp.specfit()
print "Guesses: ", sp.specfit.guesses
print "Best fit: ", sp.specfit.modelpars

raw_input('Press enter to end code')


#from matplotlib import pyplot

