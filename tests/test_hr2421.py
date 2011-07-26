import pyspeckit

if not 'interactive' in globals():
    interactive=False


sp = pyspeckit.Spectrum('hr2421.fit')

print "Does it have an axis? ",sp.plotter.axis
sp.plotter()
print "How about now? ",sp.plotter.axis

sp.plotter(xmin=4700,xmax=5000)
print "Plotter min/max: ",sp.plotter.xmin,sp.plotter.xmax," Fitter min/max: ",sp.specfit.gx1,sp.specfit.gx2," Fitregion= ",sp.baseline.fitregion


if interactive: raw_input("Wait here a moment")
sp.baseline(subtract=False,exclude=[4830,4890])
print "Plotter min/max: ",sp.plotter.xmin,sp.plotter.xmax," Fitter min/max: ",sp.specfit.gx1,sp.specfit.gx2," Fitregion= ",sp.baseline.fitregion
print "Baseline exclude: ",sp.baseline.excludevelo,sp.baseline.excludepix
if interactive: raw_input("Wait here a moment")

# set the baseline to zero to prevent variable-height fitting
# (if you don't do this, the best fit to the spectrum is dominated by the
# background level)
#sp.baseline.order = 0
print "FITTING GAUSSIAN"
sp.specfit()
sp.specfit() # Do this twice to get a better estimate of the noise
print "Plotter min/max: ",sp.plotter.xmin,sp.plotter.xmax," Fitter min/max: ",sp.specfit.gx1,sp.specfit.gx2," Fitregion= ",sp.baseline.fitregion
sp.plotter.figure.savefig('hr2421_gaussfit.png')
print "Guesses: ", sp.specfit.guesses
print "Best fit: ", sp.specfit.modelpars
print "EQW: ",sp.specfit.EQW()
print "Chi2: ",sp.specfit.chi2
gauss_model = sp.specfit.model+sp.baseline.basespec[sp.specfit.gx1:sp.specfit.gx2]
if interactive: raw_input("Wait here a moment")

print "FITTING VOIGT"
sp.specfit(fittype='voigt')
print "Guesses: ", sp.specfit.guesses
print "Best fit: ", sp.specfit.modelpars
print "EQW: ",sp.specfit.EQW()
print "Chi2: ",sp.specfit.chi2
sp.plotter.axis.plot(sp.xarr[sp.specfit.gx1:sp.specfit.gx2],gauss_model,color='b',linewidth=0.5)
sp.plotter(clear=False)
sp.plotter.figure.savefig('hr2421_voigtfit.png')
voigt_model = sp.specfit.model+sp.baseline.basespec[sp.specfit.gx1:sp.specfit.gx2]
if interactive: raw_input("Wait here a moment")

sp.specfit(interactive=True)
if interactive: raw_input('Press enter to print guesses and best fit and end code')
sp.plotter.figure.savefig('hr2421_interactive_fit.png')
print "Guesses: ", sp.specfit.guesses
print "Best fit: ", sp.specfit.modelpars

print "EQW: ",sp.specfit.EQW()


#from matplotlib import pyplot

