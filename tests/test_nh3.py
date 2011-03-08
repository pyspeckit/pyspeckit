import spectrum

sp = spectrum.Spectrum('G031.947+00.076_nh3_11_Tastar.fits')#,wcstype='F')

#sp.plotter(xmin=2.36875e10,xmax=2.36924e10) 
sp.plotter(xmin=-100,xmax=300)
#raw_input("Plotter")

# set the baseline to zero to prevent variable-height fitting
# (if you don't do this, the best fit to the spectrum is dominated by the
# background level)
sp.baseline.order = 0
sp.specfit()
sp.plotter.figure.savefig('nh3_gaussfit.png')
print "Guesses: ", sp.specfit.guesses
print "Best fit: ", sp.specfit.modelpars

#sp.baseline(exclude=[0,200],order=0)
#print "Plotter min/max: ",sp.plotter.xmin,sp.plotter.xmax," Fitter min/max: ",sp.specfit.gx1,sp.specfit.gx2," Fitregion= ",sp.baseline.excludevelo,sp.baseline.excludepix
#raw_input('Baseline')
#sp.specfit(fittype='ammonia',multifit=True,guesses=[20,20,1e16,1.0,-55.0,0.5],quiet=False,xunits='Hz')
sp.specfit(fittype='ammonia',multifit=True,guesses=[5.9,4.45,8.3e14,0.84,96.2,0.43],quiet=False,xunits=sp.xunits)
sp.specfit.plotresiduals()
raw_input('Press enter to print guesses and zoom in.')
sp.plotter.figure.savefig('nh3_ammonia_fit.png')
print "Guesses: ", sp.specfit.guesses
print "Best fit: ", sp.specfit.modelpars
sp.plotter(xmin=70,xmax=125)
sp.specfit.plot_fit()
sp.plotter.figure.savefig('nh3_ammonia_fit_zoom.png')
raw_input('Press enter to multifit')
sp.specfit(fittype='ammonia',multifit=True,
        guesses=[4,3.5,5e14,0.68,97.3,0.5]+[15,4.2,7e14,0.52,95.8,0.35],
        quiet=False,xunits=sp.xunits)
sp.specfit.plotresiduals()
sp.plotter.figure.savefig('nh3_ammonia_multifit_zoom.png')


raw_input('Press enter to end code')
