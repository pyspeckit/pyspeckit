import pyspeckit
from pylab import *

if not 'interactive' in globals():
    interactive=False
if not 'savedir' in globals():
    savedir = ''

filenames = {'oneone':'G032.751-00.071_nh3_11_Tastar.fits',
    'twotwo':'G032.751-00.071_nh3_22_Tastar.fits',
    'threethree':'G032.751-00.071_nh3_33_Tastar.fits',
    'fourfour':'G032.751-00.071_nh3_44_Tastar.fits'}

spdict1,spectra1 = pyspeckit.wrappers.fitnh3.fitnh3tkin(filenames, crop=[0,80],
        tkin=18.65, tex=4.49, column=15.5, fortho=0.9, verbose=False, smooth=6,
        fignum=6, quiet=True, rebase=True)
print spectra1.specfit.Registry
print spectra1.specfit.Registry.multifitters['ammonia']

# a sanity check
"""
line = 'oneone'
subplot(221); plot(spdict1[line].xarr,pyspeckit.models.ammonia.ammonia(spdict1[line].xarr,tkin=22.5,tex=4.5,ntot=14.5,width=1.03,xoff_v=37.94,fortho=1))
line = 'twotwo'
subplot(222); plot(spdict1[line].xarr,pyspeckit.models.ammonia.ammonia(spdict1[line].xarr,tkin=22.5,tex=4.5,ntot=14.5,width=1.03,xoff_v=37.94,fortho=1))
line = 'threethree'
subplot(223); plot(spdict1[line].xarr,pyspeckit.models.ammonia.ammonia(spdict1[line].xarr,tkin=22.5,tex=4.5,ntot=14.5,width=1.03,xoff_v=37.94,fortho=1))
line = 'fourfour'
subplot(224); plot(spdict1[line].xarr,pyspeckit.models.ammonia.ammonia(spdict1[line].xarr,tkin=22.5,tex=4.5,ntot=14.5,width=1.03,xoff_v=37.94,fortho=1))

am = pyspeckit.models.ammonia.ammonia_model()
mymodel = am.n_ammonia(pars=[22.5,4.5,14.5,1.03,37.94,1],parnames=['tkin','tex','ntot','width','xoff_v','fortho'])
line = 'oneone'
subplot(221); plot(spdict1[line].xarr,mymodel(spdict1[line].xarr))
line = 'twotwo'
subplot(222); plot(spdict1[line].xarr,mymodel(spdict1[line].xarr))
line = 'threethree'
subplot(223); plot(spdict1[line].xarr,mymodel(spdict1[line].xarr))
line = 'fourfour'
subplot(224); plot(spdict1[line].xarr,mymodel(spdict1[line].xarr))
"""

filenames_para = {'oneone':'G032.751-00.071_nh3_11_Tastar.fits',
    'twotwo':'G032.751-00.071_nh3_22_Tastar.fits',
    'fourfour':'G032.751-00.071_nh3_44_Tastar.fits'}
#spdict2,spectra2 = pyspeckit.wrappers.fitnh3.fitnh3tkin(filenames_para,crop=[0,80],tkin=18.64,tex=4.49,column=14.8,fortho=0.0,fixed=[False,False,False,False,False,True],fignum=3,guessfignum=4,smooth=6, rebase=True)
#spdict3,spectra3 = pyspeckit.wrappers.fitnh3.fitnh3tkin(filenames_para,crop=[0,80],tkin=18.64,tex=4.49,tau=4.0,fortho=0.0,fixed=[False,False,False,False,False,True],fignum=3,guessfignum=4,smooth=6,doplot=True, rebase=True)
spdict4,spectra4 = pyspeckit.wrappers.fitnh3.fitnh3tkin(filenames_para,crop=[0,80],guesses=[21.1,3.3,4.0,1.23,38.0,0.0,48.69,10.19,0.4,1.23,38.0,0.0],
        parnames=['tkin','tex','tau','width','xoff_v','fortho']*2,fixed=[False,False,False,False,False,True]*2,
        npeaks=2,fignum=3,guessfignum=4,smooth=6,doplot=True, rebase=True, tau=True, shh=False)
    # for some reason, fixed breaks multiparameter mpfit
    # fixed=[False,False,False,False,False,True]*2


if False:
    figure(7,figsize=[16,12])
    figure(8,figsize=[16,12])
    figure(5,figsize=[16,12])

    spectra1.specfit(fittype='ammonia',quiet=False,multifit=True,
            guesses=[17.57,4.36,15.49,0.82,37.96,0.86,22.49,2.97,16.06,2.19,37.88,0.84])
            #fixed=[True,False,False,True,True,False,False,False,False,False,False,False])
    spectra1.error[:] = spectra1.specfit.residuals.std()
    splist = spdict1.values()
    for sp in splist:
        sp.xarr.convert_to_unit('km/s',quiet=True)
        sp.specfit.fitter = spectra1.specfit.fitter
        sp.specfit.modelpars = spectra1.specfit.modelpars
        sp.specfit.npeaks = spectra1.specfit.npeaks
        sp.specfit.model = pyspeckit.models.ammonia.ammonia_model(npeaks=2).n_ammonia(pars=spectra1.specfit.modelpars, parnames=['tkin','tex','ntot','width','xoff_v','fortho']*2)(sp.xarr)
        sp.specfit.residuals = sp.data - sp.specfit.model
        sp.error[:] = sp.specfit.residuals.std()

    pyspeckit.wrappers.fitnh3.plot_nh3(spdict1,spectra1,fignum=5,residfignum=8)
    pyspeckit.wrappers.fitnh3.plot_nh3(spdict1,spectra1,fignum=7,show_components=True)

    if savedir != "":
        figure(5)
        savefig(savedir+"example_G032_multi-temperature_four-line_fit.png")
        figure(7)
        savefig(savedir+"example_G032_multi-temperature_four-line_fit_components.png")
        figure(8)
        savefig(savedir+"example_G032_multi-temperature_four-line_fit_residuals.png")

    figure(9,figsize=[16,12])
    pyspeckit.wrappers.fitnh3.plot_nh3(spdict1,spectra1,fignum=9,errstyle='fill')
    if savedir != "":
        savefig(savedir+"example_G032_multi-temperature_four-line_fit_errbars.png")
