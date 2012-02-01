import pyspeckit
import numpy as np
from pyspeckit.spectrum.models import voigtfitter


xarr = pyspeckit.spectrum.units.SpectroscopicAxis(np.linspace(-100,100,500),unit='km/s',refX=1e9,refX_units='Hz')
VF = voigtfitter.voigt_fitter()

sp = pyspeckit.Spectrum(xarr=xarr, data=VF.voigt(xarr,1,0,2.5,2.5) + np.random.randn(xarr.shape[0])/20., error=np.ones(xarr.shape[0])/20.)
sp.plotter()
sp.specfit(fittype='voigt')
sp.specfit(fittype='lorentzian',composite_fit_color='g',clear=False,annotate=False)
sp.specfit(fittype='gaussian',composite_fit_color='b',clear=False,annotate=False)

sp = pyspeckit.Spectrum(xarr=xarr, data=VF.voigt(xarr,1,0,2.5,5.0) + np.random.randn(xarr.shape[0])/20., error=np.ones(xarr.shape[0])/20.)
sp.plotter()
sp.specfit(fittype='voigt')
sp.specfit(fittype='lorentzian',composite_fit_color='g',clear=False,annotate=False)
sp.specfit(fittype='gaussian',composite_fit_color='b',clear=False,annotate=False)

sp = pyspeckit.Spectrum(xarr=xarr, data=VF.voigt(xarr,1,0,2.5,5.0) + np.random.randn(xarr.shape[0])/50., error=np.ones(xarr.shape[0])/50.)
sp.plotter()
sp.specfit(fittype='voigt')
sp.specfit(fittype='lorentzian',composite_fit_color='g',clear=False,annotate=False)
sp.specfit(fittype='gaussian',composite_fit_color='b',clear=False,annotate=False)
