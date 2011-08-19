import pyspeckit
from pylab import *

filenames = {'oneone':'G032.751-00.071_nh3_11_Tastar.fits',
    'twotwo':'G032.751-00.071_nh3_22_Tastar.fits',
    'threethree':'G032.751-00.071_nh3_33_Tastar.fits',
    'fourfour':'G032.751-00.071_nh3_44_Tastar.fits'}

spdict1,spectra1 = pyspeckit.wrappers.fitnh3.fitnh3tkin(filenames,crop=[0,80],tkin=18.65,tex=4.49,column=14.8,fortho=0.7,verbose=False,smooth=5, fignum=6)

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
spdict2,spectra2 = pyspeckit.wrappers.fitnh3.fitnh3tkin(filenames_para,crop=[0,80],tkin=18.64,tex=4.49,column=14.8,fortho=0.7,fignum=3,guessfignum=4,smooth=5)



spectra1.specfit(fittype='ammonia',quiet=False,multifit=True,
        guesses=[18.64,4.49,14.5,1.03,37.94,0.7,30.3,5.1,15.0,4.0,38.0,0.3],
        fixed=[True,False,False,True,True,False,False,False,False,False,False,False])
splist = spdict1.values()
for sp in splist:
    sp.xarr.convert_to_unit('km/s',quiet=True)
    sp.specfit.fitter = spectra1.specfit.fitter
    sp.specfit.modelpars = spectra1.specfit.modelpars
    sp.specfit.npeaks = spectra1.specfit.npeaks
    sp.specfit.model = pyspeckit.models.ammonia.ammonia_model(npeaks=2).n_ammonia(pars=spectra1.specfit.modelpars, parnames=['tkin','tex','ntot','width','xoff_v','fortho']*2)(sp.xarr)

pyspeckit.wrappers.fitnh3.plot_nh3(spdict1,spectra1,fignum=5)
