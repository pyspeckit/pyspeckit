import pyspeckit
import numpy as np
from pyspeckit.spectrum.models import inherited_voigtfitter

# This example uses scipy
try:
    import scipy
except ImportError:
    exit

# technically, the voigt fitter works as a singlefitter (i.e., you can fit the
# background level and the peak simultaneously)
# in practice, however, you need to fit the background independently except for
# gaussians.  I don't know why this is.

xarr = pyspeckit.spectrum.units.SpectroscopicAxis(np.linspace(-100, 100, 500),
                                                  unit='km/s',
                                                  refX=1e9,
                                                  refX_unit='Hz')
VF = inherited_voigtfitter.voigt_fitter()

sp1 = pyspeckit.Spectrum(xarr=xarr,
                         data=(VF.n_modelfunc((1, 0, 2.5, 2.5))(xarr) +
                               np.random.randn(xarr.shape[0])/20.),
                         error=np.ones(xarr.shape[0])/20.,
                         header={},
                        )
sp1.plotter()
sp1.specfit(fittype='gaussian', composite_fit_color='b', clear=False,
            annotate=False, guesses='moments')
sp1.specfit(fittype='lorentzian', composite_fit_color='g', clear=False,
            annotate=False, guesses='moments')
sp1.specfit(fittype='voigt', composite_fit_color='r', clear=False,
            annotate=True, guesses='moments')

sp2 = pyspeckit.Spectrum(xarr=xarr, data=VF.n_modelfunc((1,0,2.5,5.0))(xarr) +
                         np.random.randn(xarr.shape[0])/20.,
                         error=np.ones(xarr.shape[0])/20.,
                         header={},
                        )
sp2.plotter()
sp2.specfit(fittype='gaussian', composite_fit_color='b', clear=False,
            annotate=False, guesses='moments')
sp2.specfit(fittype='lorentzian', composite_fit_color='g', clear=False,
            annotate=False, guesses='moments')
sp2.specfit(fittype='voigt', composite_fit_color='r', clear=False,
            annotate=True, guesses='moments')

sp3 = pyspeckit.Spectrum(xarr=xarr, data=VF.n_modelfunc((1,0,2.5,5.0))(xarr) +
                         np.random.randn(xarr.shape[0])/50.,
                         error=np.ones(xarr.shape[0])/50.,
                         header={},
                        )
sp3.plotter()
sp3.specfit(fittype='gaussian', composite_fit_color='b', clear=False,
            annotate=False, guesses='moments')
sp3.specfit(fittype='lorentzian', composite_fit_color='g', clear=False,
            annotate=False, guesses='moments')
sp3.specfit(fittype='voigt', composite_fit_color='r', clear=False,
            annotate=True, guesses='moments')
