import numpy as np
import warnings

from pyspeckit.spectrum import Spectrum


def test_eqw():
    dx = 0.1
    x = np.arange(-6,6,dx)
    y = 1-np.exp(-x**2 / 2.)
    with warnings.catch_warnings():
        # ignore warning about creating an empty header
        warnings.simplefilter('ignore')
        sp = Spectrum(xarr=x, data=y)
    sp.baseline(exclude=[-5,5], order=0, subtract=False)
    sp.specfit(fittype='gaussian', guesses=(-1, 0, 0.5))
    eqw_nofit = sp.specfit.EQW(fitted=False)
    np.testing.assert_almost_equal(sp.specfit.EQW(), (1-y).sum() * dx, decimal=5)
    np.testing.assert_almost_equal(sp.specfit.EQW(continuum=1), (1-y).sum() * dx, decimal=5)
    np.testing.assert_almost_equal(eqw_nofit, (1-y).sum() * dx, decimal=4)
    return sp

def test_eqw_plot():
    dx = 0.1
    x = np.arange(-6,6,dx)
    y = 1-np.exp(-x**2 / 2.)
    with warnings.catch_warnings():
        # ignore warning about creating an empty header
        warnings.simplefilter('ignore')
        sp = Spectrum(xarr=x, data=y)
    sp.plotter()
    sp.baseline(exclude=[-5,5], order=0, subtract=False)
    sp.specfit(fittype='gaussian', guesses=(-1, 0, 0.5))
    eqw = sp.specfit.EQW()
    eqw_nofit = sp.specfit.EQW(fitted=False)
    eqw_cont = sp.specfit.EQW(plotcolor='g', fitted=True, continuum=1,
                              plot=True,
                              components=False, annotate=True,
                              loc='lower left', xmin=None, xmax=None)
    np.testing.assert_almost_equal(eqw, (1-y).sum() * dx, decimal=5)
    np.testing.assert_almost_equal(eqw_cont, (1-y).sum() * dx, decimal=5)
    np.testing.assert_almost_equal(eqw_nofit, (1-y).sum() * dx, decimal=4)
    return sp
