from .. import hill5infall
from ...classes import Spectrum
from ... import units
import numpy as np
from astropy import units as u

def test_hill5():
    x = np.linspace(-5,5,20)*u.km/u.s
    xarr = units.SpectroscopicAxis(x, refX=5*u.GHz, velocity_convention='radio')
    y = hill5infall.hill5_model(xarr, 0.5, 1.0, 2.0, 1.0, 1.0)

    sp = Spectrum(xarr=xarr, data=y, header={})

    sp.Registry.add_fitter('hill5_fitter', hill5infall.hill5_fitter, 5)

    sp.specfit(fittype='hill5_fitter', guesses=[0.4, 0.9, 2.2, 0.9, 1.1])

    np.testing.assert_almost_equal(sp.specfit.parinfo[0].value, 0.5)
    np.testing.assert_almost_equal(sp.specfit.parinfo[1].value, 1.0)
    np.testing.assert_almost_equal(sp.specfit.parinfo[2].value, 2.0)
    np.testing.assert_almost_equal(sp.specfit.parinfo[3].value, 1.0)
    np.testing.assert_almost_equal(sp.specfit.parinfo[4].value, 1.0)

def test_hill5_guessing():

    rslt = hill5infall.parse_3par_guesses(['amp', 'cen', 'wid'])

    assert rslt == [1.0, 'cen', 1.0, 'wid', 'amp']

    rslt = hill5infall.parse_3par_guesses(['amp1', 'cen1', 'wid1',
                                           'amp2', 'cen2', 'wid2',
                                          ])

    assert rslt == [1.0, 'cen1', 1.0, 'wid1', 'amp1',
                    1.0, 'cen2', 1.0, 'wid2', 'amp2',
                   ]
