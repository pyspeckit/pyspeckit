import spectrum

sp = spectrum.Spectrum('sample_13CO.fits')

print "Does it have an axis? ",sp.plotter.axis
sp.plotter()
print "How about now? ",sp.plotter.axis

raw_input('blah')

#from matplotlib import pyplot

