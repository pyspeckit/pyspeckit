import pyspeckit
import matplotlib

if not 'interactive' in globals():
    interactive=False
if not 'savedir' in globals():
    savedir = ''


sp = pyspeckit.Spectrum('hr2421.fit')

print "Does it have an axis? ",sp.plotter.axis
sp.plotter()
print "How about now? ",sp.plotter.axis

sp.plotter(xmin=4700,xmax=5000)
print "Plotter min/max: ",sp.plotter.xmin,sp.plotter.xmax," Fitter min/max: ",sp.specfit.xmin,sp.specfit.xmax," Fitregion= ",sp.baseline.button1plot," bfit target sum: ",sp.baseline.includemask.sum()

if interactive: raw_input("Wait here a moment")
sp.baseline(subtract=False,exclude=[4830,4890],order=2,highlight=True)
print "Plotter min/max: ",sp.plotter.xmin,sp.plotter.xmax," Fitter min/max: ",sp.specfit.xmin,sp.specfit.xmax," Fitregion= ",sp.baseline.button1plot," bfit target sum: ",sp.baseline.includemask.sum()
# obsolete 1/16/2012 print "Baseline exclude: ",sp.baseline.excludevelo,sp.baseline.excludepix
if interactive: raw_input("Wait here a moment")

# set the baseline to zero to prevent variable-height fitting
# (if you don't do this, the best fit to the spectrum is dominated by the
# background level)
#sp.baseline.order = 0
print "FITTING GAUSSIAN"
sp.specfit(debug=True,verbose=True)
sp.specfit(debug=True,verbose=True) # Do this twice to get a better estimate of the noise
print "Plotter min/max: ",sp.plotter.xmin,sp.plotter.xmax," Fitter min/max: ",sp.specfit.xmin,sp.specfit.xmax," Fitregion= ",sp.baseline.button1plot," bfit target sum: ",sp.baseline.includemask.sum()
sp.plotter.figure.savefig(savedir+'hr2421_gaussfit.png')
print "Guesses: ", sp.specfit.guesses
print "Best fit: ", sp.specfit.modelpars
print "EQW: ",sp.specfit.EQW()
print "Chi2: ",sp.specfit.chi2
gauss_model = sp.specfit.model+sp.baseline.basespec[sp.specfit.xmin:sp.specfit.xmax]
print "A Gaussian has been fit.  There should be a red line overlaid on the spectrum with green highlighting the baseline fit region and yellow showing the baseline fit"
if interactive: raw_input("Wait here a moment")

try:
    print "FITTING VOIGT"
    sp.specfit(fittype='voigt')
    print "Guesses: ", sp.specfit.guesses
    print "Best fit: ", sp.specfit.modelpars
    print "EQW: ",sp.specfit.EQW()
    print "Chi2: ",sp.specfit.chi2
    sp.plotter.axis.plot(sp.xarr[sp.specfit.xmin:sp.specfit.xmax],gauss_model,color='b',linewidth=0.5)
    sp.plotter(clear=False)
    sp.plotter.figure.savefig(savedir+'hr2421_voigtfit.png')
    voigt_model = sp.specfit.model+sp.baseline.basespec[sp.specfit.xmin:sp.specfit.xmax]
    print "A voigt model has been fit.  The red line from before should have a blue line overlaid.  They should be only moderately different."
    if interactive: raw_input("Wait here a moment")
except ImportError:
    print "Could not fit voigt profiles because scipy wasn't installed."


sp.specfit(interactive=True)
if interactive:
    raw_input('Press enter to print guesses and best fit and end code')
else:
    event1 = matplotlib.backend_bases.MouseEvent('button_press_event', sp.plotter.axis.figure.canvas,257,316,button=1)
    event2 = matplotlib.backend_bases.MouseEvent('button_press_event', sp.plotter.axis.figure.canvas,732,280,button=1)
    event3 = matplotlib.backend_bases.MouseEvent('button_press_event', sp.plotter.axis.figure.canvas,523,194,button=2)
    event4 = matplotlib.backend_bases.MouseEvent('button_press_event', sp.plotter.axis.figure.canvas,485,264,button=2)
    event5 = matplotlib.backend_bases.MouseEvent('button_press_event', sp.plotter.axis.figure.canvas,611,247,button=3)
    sp.specfit.event_manager(event1)
    sp.specfit.event_manager(event2)
    sp.plotter.figure.savefig(savedir+'hr2421_interactive_selectregion.png')
    sp.specfit.event_manager(event3)
    sp.specfit.event_manager(event4)
    sp.plotter.figure.savefig(savedir+'hr2421_interactive_guesses.png')
    sp.specfit.event_manager(event5)

sp.plotter.figure.savefig(savedir+'hr2421_interactive_fit.png')
print "Guesses: ", sp.specfit.guesses
print "Best fit: ", sp.specfit.modelpars

print "EQW: ",sp.specfit.EQW()


#from matplotlib import pyplot

