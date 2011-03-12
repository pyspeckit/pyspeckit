import spectrum

sp1 = spectrum.Spectrum('G032.751-00.071_nh3_11_Tastar.fits')
sp1.crop(0,70)
sp2 = spectrum.Spectrum('G032.751-00.071_nh3_22_Tastar.fits')
sp2.crop(0,70)
#sp4 = spectrum.Spectrum('G032.751-00.071_nh3_44_Tastar.fits')
#sp4.crop(3460,3500)
spectra = spectrum.Spectra([sp1,sp2])
sp = spectra

sp.plotter()
#sp.plotter(xmin=2.36875e10,xmax=2.36924e10) 
#sp.plotter(xmin=-100,xmax=300)
#raw_input("Plotter")
from pylab import *
show()


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
sp.specfit(fittype='ammonia', multifit=True, guesses=[21, 0.3, 5.3e13, 1.14,
    35, 0.5], fixed=[False,False,False,False,False,True],
    minpars=[2.73,0,1e10,0.1,0,0],
    limitedmin=[True,True,True,True,False,True], quiet=False,
    xunits=sp.xarr.units)
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
        guesses=[25,0.5,5e14,0.68,37,0.5]+[25,0.2,7e14,0.52,37,0.5],
        fixed=[False,False,False,False,False,True]*2,
        maxpars=[50,0,0,0,0,1]*2,
        limitedmax=[True,False,False,False,False,True]*2,
        minpars=[2.73,0,1e10,0.1,0,0]*2,
        limitedmin=[True,True,True,True,False,True]*2,
        quiet=False,xunits=sp.xarr.units)
sp.specfit.plotresiduals()
sp.plotter.figure.savefig('nh3_ammonia_multifit_zoom.png')



raw_input('Press enter to end code')
