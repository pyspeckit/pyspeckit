from pyspeckit.spectrum.readers import gbt
session = gbt.GBTSession('3C286.fits')
print session
session.load_target('3C286')
target = session['3C286']
sp = target.blocks['A13OFF2'][0] # there are 4 identical spectra
sp.xarr.convert_to_unit('km/s')
sp.plotter(xmin=207408,xmax=207598)
sp.plotter.label(verbose_label=True)
stats = sp.stats(statrange=(-20+207458,20+207458))
sp.error[:] = stats['std']
sp.specfit(fittype='gaussian')
sp.specfit.plotresiduals(drawstyle='line')
print "Gaussian chi^2: %g  chi^2/n: %g" % (sp.specfit.chi2, sp.specfit.chi2/sp.specfit.dof)
print "Optimal chi^2: %g  chi^2/n: %g" % (sp.specfit.optimal_chi2(reduced=False),sp.specfit.optimal_chi2())
sp.specfit(fittype='voigt', clear=False, composite_fit_color='blue')
sp.specfit.plotresiduals(clear=False, color='blue', drawstyle='line')
print "Voigt   chi^2: %g  chi^2/n: %g" % (sp.specfit.chi2, sp.specfit.chi2/sp.specfit.dof)
print "Optimal chi^2: %g  chi^2/n: %g" % (sp.specfit.optimal_chi2(reduced=False),sp.specfit.optimal_chi2())
sp.specfit.residualaxis.set_ylim(-0.2,0.2)

# this is usually unnecessary, but it's platform-dependent
import pylab
pylab.ion()
pylab.show()


import pyspeckit
sp1 = target.blocks['A13OFF2'][0]
sp2 = target.blocks['A13OFF2'][1]
xc = pyspeckit.correlate.correlate(sp1,sp2,range=[207408,207698],units='km/s')
xc.plotter()
xc.specfit()

