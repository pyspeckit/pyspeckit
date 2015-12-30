from .. import hill5infall
from ...classes import Spectrum
from ... import units
import numpy as np
from astropy import units as u

def test_hill5():
    x = np.linspace(-5,5,20)*u.km/u.s
    xarr = units.SpectroscopicAxis(x, refX=5*u.GHz, velocity_convention='radio')
    y = hill5infall.hill5_model(xarr, 0.5, 1.0, 2.0, 1.0, 1.0)

    sp = Spectrum(xarr=xarr, data=y)

    sp.Registry.add_fitter('hill5_fitter', hill5infall.hill5_fitter, 5)

    sp.specfit(fittype='hill5_fitter', guesses=[0.4, 0.9, 2.2, 0.9, 1.1])

    np.testing.assert_almost_equal(sp.specfit.parinfo[0].value, 0.5)
    np.testing.assert_almost_equal(sp.specfit.parinfo[1].value, 1.0)
    np.testing.assert_almost_equal(sp.specfit.parinfo[2].value, 2.0)
    np.testing.assert_almost_equal(sp.specfit.parinfo[3].value, 1.0)
    np.testing.assert_almost_equal(sp.specfit.parinfo[4].value, 1.0)
