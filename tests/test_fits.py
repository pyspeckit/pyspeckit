import spectrum

sp = spectrum.Spectrum('sample_13CO.fits')

print "Does it have an axis? ",sp.plotter.axis
sp.plotter()
print "How about now? ",sp.plotter.axis


sp.baseline.order = 0
sp.specfit(quiet=False)
print "Guesses: ", sp.specfit.guesses
print "Best fit: ", sp.specfit.modelpars

raw_input('Press enter to end code')


#from matplotlib import pyplot

