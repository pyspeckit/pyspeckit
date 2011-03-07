import spectrum

sp = spectrum.Spectrum('simple_txt.txt')

print "Does it have an axis? ",sp.plotter.axis
sp.plotter()
print "How about now? ",sp.plotter.axis


print "FITTING GAUSSIAN"
sp.specfit(quiet=False)
sp.specfit.annotate(loc='lower right')
sp.plotter.figure.savefig('txt_gaussfit.png')
print "Guesses: ", sp.specfit.guesses
print "Best fit: ", sp.specfit.modelpars
print "Datamax: ",sp.data.max()
raw_input("Wait here a moment")

# Fit a baseline, excluding the velocities with data, and don't subtract it
print "FITTING BASELINE"
sp.baseline(exclude=[2.111,2.129],subtract=False)
sp.plotter.figure.savefig('txt_baseline.png')
print "Baseline: ",sp.baseline.baselinepars
print "Excludepix: ",sp.baseline.excludepix
print "Excludevelo: ",sp.baseline.excludevelo
print "EQW: ",sp.specfit.EQW()
print "Datamax: ",sp.data.max()
print "NOK: ",sp.baseline.OKmask.sum()
print sp.data[True-sp.baseline.excludemask]
print sp.baseline.basespec
raw_input("Wait here a moment")

print "REFITTING GAUSSIAN"
sp.specfit(quiet=False)
sp.specfit.annotate(loc='lower right')
print "Guesses: ", sp.specfit.guesses
print "Best fit: ", sp.specfit.modelpars
print "Datamax: ",sp.data.max()
sp.plotter.figure.savefig('txt_baseline_gaussfit.png')
raw_input("Wait here a moment")

print "EQW: ",sp.specfit.EQW(plot=True,annotate=True)
sp.plotter.refresh()
sp.plotter.figure.savefig('txt_EQW.png')

raw_input('Press enter to print guesses and best fit and end code')


#from matplotlib import pyplot


