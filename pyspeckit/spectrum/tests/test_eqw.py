import numpy as np
import warnings

import pytest

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

def test_eqw_fitted_midpt():
    # Regression test for issue #96: EQW(plot=True, midpt_location='fitted')
    # used to crash (undefined variable `sp` and unset midpt_pixel)
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
    eqw_default = sp.specfit.EQW()
    eqw_fitted = sp.specfit.EQW(plot=True, midpt_location='fitted')
    assert np.isfinite(eqw_fitted)
    # midpt_location only affects where the EQW box is drawn, not the value
    np.testing.assert_almost_equal(eqw_fitted, eqw_default, decimal=6)
    np.testing.assert_almost_equal(eqw_fitted, (1-y).sum() * dx, decimal=5)
    # the plotted EQW rectangle should be centered on the fitted line center
    verts = sp.specfit.EQW_plots[-1].get_paths()[0].vertices
    box_center = (verts[:,0].min() + verts[:,0].max()) / 2.
    np.testing.assert_almost_equal(box_center,
                                   sp.specfit.parinfo['SHIFT0'].value,
                                   decimal=6)
    return sp

def test_eqw_fitted_midpt_no_shift_param():
    # Regression test for issue #96: a model without any SHIFT parameter
    # should raise an informative AttributeError, not an obscure ValueError
    from pyspeckit.spectrum.models import model

    def offset_gaussian(x, amp, xoff_v, width):
        return amp * np.exp(-(x - xoff_v)**2 / (2. * width**2))

    fitter = model.SpectralModel(offset_gaussian, 3,
                                 parnames=['amplitude', 'xoff_v', 'width'],
                                 shortvarnames=('A', 'v', 'w'))
    dx = 0.1
    x = np.arange(-6,6,dx)
    y = 1-np.exp(-x**2 / 2.)
    with warnings.catch_warnings():
        # ignore warning about creating an empty header
        warnings.simplefilter('ignore')
        sp = Spectrum(xarr=x, data=y)
    sp.plotter()
    sp.baseline(exclude=[-5,5], order=0, subtract=False)
    sp.specfit.register_fitter('offset_gaussian', fitter, 3)
    sp.specfit(fittype='offset_gaussian', guesses=(-1, 0, 0.5))
    with pytest.raises(AttributeError):
        sp.specfit.EQW(plot=True, midpt_location='fitted')
