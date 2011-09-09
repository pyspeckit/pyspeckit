import pyspeckit
import numpy as np


sp1 = pyspeckit.Spectrum('G203.04+1.76_h2co.fits',wcstype='D',scale_keyword='ETAMB')
#sp1.data /= sp1.header.get('ETAMB')
sp2 = pyspeckit.Spectrum('G203.04+1.76_h2co_Tastar.fits',wcstype='V',scale_keyword='ETAMB')
#sp2.data /= sp2.header.get('ETAMB')

# pyspeckit.classes.register_fitter(sp1.Registry,'formaldehyde_radex',
#         pyspeckit.models.formaldehyde.formaldehyde11_radex_fitter,4,multisingle='multi')
# 
# sp1.plotter()
# sp1.specfit(fittype='formaldehyde_radex',multifit=True,guesses=[4,12,3.75,0.43],quiet=False)

sp1.crop(-50,50)
sp1.smooth(3) # match to GBT resolution
sp2.crop(-50,50)

sp1.xarr.convert_to_unit('GHz')
sp1.specfit() # determine errors
sp1.error = np.ones(sp1.data.shape)*sp1.specfit.residuals.std()
sp2.xarr.convert_to_unit('GHz')
sp2.specfit() # determine errors
sp2.error = np.ones(sp2.data.shape)*sp2.specfit.residuals.std()
sp = pyspeckit.Spectra([sp1,sp2])

pyspeckit.classes.register_fitter(sp.Registry,'formaldehyde_radex_2',
        pyspeckit.models.formaldehyde.formaldehyde_radex_fitter,4,multisingle='multi')

sp.plotter()
sp.specfit(fittype='formaldehyde_radex_2',multifit=True,guesses=[4,12,3.75,0.43],quiet=False)

sp1.specfit.fitter = sp.specfit.fitter
sp1.specfit.model = np.interp(sp1.xarr,sp.xarr,sp.specfit.model)
sp2.specfit.fitter = sp.specfit.fitter
sp2.specfit.model = np.interp(sp2.xarr,sp.xarr,sp.specfit.model)

sp1.xarr.convert_to_unit('km/s')
sp2.xarr.convert_to_unit('km/s')

sp1.plotter(xmin=-5,xmax=15,errstyle='fill')
sp1.specfit.plot_fit()
sp2.plotter(xmin=-5,xmax=15,errstyle='fill')
sp2.specfit.plot_fit()
