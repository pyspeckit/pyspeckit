import pyspeckit

if not 'interactive' in globals():
    interactive=False
if not 'savedir' in globals():
    savedir = ''

filenames = {'oneone':'G032.751-00.071_nh3_11_Tastar.fits',
    'twotwo':'G032.751-00.071_nh3_22_Tastar.fits',
    'threethree':'G032.751-00.071_nh3_33_Tastar.fits',}

sp1 = pyspeckit.Spectrum('G032.751-00.071_nh3_11_Tastar.fits')
sp1.crop(0,80)
sp1.smooth(4)
sp1.plotter()
sp1.baseline(exclude=[15,62])
sp1.plotter()
sp1.plotter.savefig(savedir+'nh3_11_baselined.png')
sp2 = pyspeckit.Spectrum('G032.751-00.071_nh3_22_Tastar.fits')
sp2.crop(0,80)
sp2.smooth(4)
sp2.plotter()
sp2.specfit()
sp2.baseline(exclude=[34,42])
sp2.plotter()
sp2.plotter.savefig(savedir+'nh3_22_baselined.png')
sp3 = pyspeckit.Spectrum('G032.751-00.071_nh3_33_Tastar.fits')
sp3.crop(0,80)
sp3.smooth(4)
sp3.plotter()
sp3.specfit()
sp3.baseline(exclude=[30,42])
sp3.plotter()
sp3.plotter.savefig(savedir+'nh3_33_baselined.png')
#sp4.crop(2.3868e10,2.3871e10)
spectra = pyspeckit.Spectra([sp1,sp2,sp3])
sp = spectra

sp.plotter()
#sp.plotter(xmin=2.36875e10,xmax=2.36924e10) 
#sp.plotter(xmin=-100,xmax=300)
#if interactive: raw_input("Plotter")
from pylab import *
draw()
#if interactive: raw_input('wait for plotter')

ammonia_model = pyspeckit.models.ammonia_model()
published_model = pyspeckit.models.ammonia.ammonia(sp1.xarr,tkin=21.57,tex=0.24+2.73,width=1.11,xoff_v=37.88,tau=3.06)

# set the baseline to zero to prevent variable-height fitting
# (if you don't do this, the best fit to the spectrum is dominated by the
# background level)
sp.baseline.order = 0
sp.specfit()
sp.plotter.figure.savefig(savedir+'nh3_gaussfit.png')
print "Guesses: ", sp.specfit.guesses
print "Best fit: ", sp.specfit.modelpars

#sp.baseline(exclude=[0,200],order=0)
#print "Plotter min/max: ",sp.plotter.xmin,sp.plotter.xmax," Fitter min/max: ",sp.specfit.gx1,sp.specfit.gx2," Fitregion= ",sp.baseline.excludevelo,sp.baseline.excludepix
#if interactive: raw_input('Baseline')
#sp.specfit(fittype='ammonia',multifit=True,guesses=[20,20,1e16,1.0,-55.0,0.5],quiet=False,xunits='Hz')
sp.specfit(fittype='ammonia', multifit=True, guesses=[21.57, 5.0, 14.469, 1.11,
    37.8, 0.5], fixed=[False,False,False,False,False,True],
    minpars=[2.73,2.73,10,0.1,0,0],
    limitedmin=[True,True,True,True,False,True], quiet=False,
    )
sp.specfit.plotresiduals()
sp.plotter.figure.savefig(savedir+'nh3_ammonia_multifit.png')
print "Guesses: ", sp.specfit.guesses
print "Best fit: ", sp.specfit.modelpars
sp.specfit.plot_fit()
sp.plotter(xmin=2.3689e10,xmax=2.3694e10,clear=False)
sp.plotter.figure.savefig(savedir+'nh3_ammonia_multifit_11_zoom.png')
sp.plotter(xmin=2.3717e10,xmax=2.3722e10,clear=False)
sp.plotter.figure.savefig(savedir+'nh3_ammonia_multifit_22_zoom.png')
sp.plotter(xmin=2.387013e10-5.4e6,xmax=2.387013e10,clear=False)
sp.plotter.figure.savefig(savedir+'nh3_ammonia_multifit_33_zoom.png')
sp.xarr.convert_to_unit('GHz',frame='LSR')
sp.plotter(axis=subplot(131),xmin=2.3689e1,xmax=2.3694e1)
sp.specfit.plot_fit()
sp.plotter(axis=subplot(132),xmin=2.3717e1,xmax=2.3722e1)
sp.specfit.plot_fit()
sp.plotter.axis.yaxis.label.set_visible(False)
sp.plotter(axis=subplot(133),xmin=2.387013e1-5.4e-3,xmax=2.387013e1)
sp.specfit.plot_fit()
sp.plotter.axis.yaxis.label.set_visible(False)
sp.plotter.figure.savefig(savedir+'nh3_ammonia_multifit_multipanel_zoom.png')
"""
if interactive: raw_input('Press enter to multifit')
sp.specfit(fittype='ammonia',multifit=True,
        guesses=[25,0.5,5e14,0.68,37,0.5]+[25,0.2,7e14,0.52,37,0.5],
        fixed=[False,False,False,False,False,True]*2,
        maxpars=[50,0,0,0,0,1]*2,
        limitedmax=[True,False,False,False,False,True]*2,
        minpars=[2.73,0,10,0.1,0,0]*2,
        limitedmin=[True,True,True,True,False,True]*2,
        quiet=False,)
sp.specfit.plotresiduals()
sp.plotter.figure.savefig(savedir+'nh3_ammonia_multifit_zoom.png')
"""

sp3.plotter()
sp3.baseline()
sp3.specfit(fittype='ammonia', multifit=True, guesses=[145.0, 135.0, 14.74, 1.752,
    36.04, 0.5], fixed=[False,False,False,False,False,True],
    minpars=[2.73,2.73,10,0.1,0,0],
    limitedmin=[True,True,True,True,False,True], quiet=False,
    )
sp3.plotter.savefig(savedir+'nh3_33_test1.png')
sp3.baseline(excludefit=True)
sp3.plotter()
if sp.specfit.modelpars[1] > 2.73+5 and sp[1].specfit.modelpars[1] < sp.specfit.modelpars[0]: sp3.specfit.modelpars[1] += 5
sp3.specfit(fittype='ammonia', multifit=True, guesses=sp3.specfit.modelpars,
    fixed=[False,False,False,False,False,True],
    minpars=[2.73,2.73,10,0.1,0,0],
    limitedmin=[True,True,True,True,False,True], quiet=False,
    )
sp3.plotter.savefig(savedir+'nh3_33_test2_rebased.png')
if sp3.specfit.modelpars[1] > 2.75+5: sp3.specfit.modelpars[1] -= 5
sp3.specfit(fittype='ammonia', multifit=True, guesses=sp3.specfit.modelpars,
    fixed=[False,False,False,False,False,False],
    minpars=[2.73,2.73,10,0.1,0,0],
    limitedmin=[True,True,True,True,False,True], quiet=False,
    )
sp3.plotter.savefig(savedir+'nh3_33_test3_freeorthopara.png')
sp3.specfit.plotresiduals(figure=5)
sp3.specfit.residualaxis.figure.savefig(savedir+'nh3_33_test3_freeorthopara_residuals.png')


sp.plotter(reset_xlimits=True,reset_ylimits=True)
sp.specfit(fittype='ammonia', multifit=True, 
    guesses=[18.2, 17.5, 15.08, 1.01, 37.95, 0.52, 32, 30, 14.41, 2.76, 38.08, 0.45],
    fixed=[False,False,False,False,False,False]*2,
    minpars=[2.73,2.73,10,0.1,0,0]+[2.73,2.73,10,0.1,0,0],
    limitedmin=[True,True,True,True,False,True]*2, quiet=False,
    show_components=True
    )
if True:
    figure(1); clf();
    sp.specfit.fullsizemodel()
    sp.xarr.convert_to_unit('GHz',frame='LSR')
    sp.plotter(axis=subplot(211),xmin=2.3689e1,xmax=2.3694e1)
    sp.specfit.plot_fit(show_components=True)
    sp.specfit.fullsizemodel()
    sp.plotter(axis=subplot(223),xmin=2.3717e1,xmax=2.3722e1)
    sp.specfit.plot_fit(show_components=True)
    sp.specfit.fullsizemodel()
    sp.plotter(axis=subplot(224),xmin=2.387013e1-5.4e-3,xmax=2.387013e1)
    sp.specfit.plot_fit(show_components=True)
    sp.specfit.annotate(labelspacing=0.05,prop={'size':'small','stretch':'extra-condensed'},frameon=False,)
    sp.plotter.figure.savefig(savedir+'nh3_ammonia_multifit_multipanel_zoom_basedon33.png')



if interactive: raw_input('Press enter to end code')
