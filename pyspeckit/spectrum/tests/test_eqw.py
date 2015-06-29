import numpy as np

from .. import Spectrum


def test_eqw():
    dx = 0.1
    x = np.arange(-6,6,dx)
    y = 1-np.exp(-x**2 / 2.)
    sp = Spectrum(xarr=x, data=y)
    sp.baseline(exclude=[-5,5], order=0, subtract=False)
    sp.specfit(fittype='gaussian', guesses=(-1, 0, 0.5))
    np.testing.assert_almost_equal(sp.specfit.EQW(), (1-y).sum() * dx, decimal=5)
    np.testing.assert_almost_equal(sp.specfit.EQW(continuum=1), (1-y).sum() * dx, decimal=5)
