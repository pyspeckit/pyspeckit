"""
Regression tests for issue #401: profile functions documented as accepting
numpy.ndarray inputs must actually accept them (not just SpectroscopicAxis /
Quantity inputs).
"""
import numpy as np
import pytest
from astropy import units as u

from pyspeckit.spectrum.units import SpectroscopicAxis
from pyspeckit.spectrum.models.inherited_voigtfitter import voigt
from pyspeckit.spectrum.models.inherited_gaussfitter import gaussian
from pyspeckit.spectrum.models.inherited_lorentzian import lorentzian

XVALS = np.linspace(-5, 5, 101)


def test_voigt_accepts_plain_ndarray():
    # regression test for #401: this used to raise
    # AttributeError: 'numpy.ndarray' object has no attribute 'value'
    result = voigt(XVALS, 1.0, 0.0, 1.0, 1.0)
    assert isinstance(result, np.ndarray)
    assert np.all(np.isfinite(result))


@pytest.mark.parametrize('normalized', (False, True))
def test_voigt_ndarray_matches_spectroscopicaxis(normalized):
    xarr = SpectroscopicAxis(XVALS, unit=u.km/u.s)
    from_axis = voigt(xarr, 1.0, 0.0, 1.0, 1.0, normalized=normalized)
    from_ndarray = voigt(XVALS, 1.0, 0.0, 1.0, 1.0, normalized=normalized)
    assert np.all(np.isfinite(from_ndarray))
    np.testing.assert_allclose(from_ndarray, from_axis)


def test_voigt_ndarray_matches_quantity():
    from_quantity = voigt(XVALS*u.GHz, 1.0, 0.0, 1.0, 1.0)
    from_ndarray = voigt(XVALS, 1.0, 0.0, 1.0, 1.0)
    np.testing.assert_allclose(from_ndarray, from_quantity)


def test_gaussian_ndarray_matches_spectroscopicaxis():
    xarr = SpectroscopicAxis(XVALS, unit=u.km/u.s)
    from_axis = gaussian(xarr, 1.0, 0.0, 1.0)
    from_ndarray = gaussian(XVALS, 1.0, 0.0, 1.0)
    assert np.all(np.isfinite(from_ndarray))
    np.testing.assert_allclose(from_ndarray, from_axis)


def test_lorentzian_ndarray_matches_spectroscopicaxis():
    xarr = SpectroscopicAxis(XVALS, unit=u.km/u.s)
    from_axis = lorentzian(xarr, 1.0, 0.0, 1.0)
    from_ndarray = lorentzian(XVALS, 1.0, 0.0, 1.0)
    assert np.all(np.isfinite(from_ndarray))
    np.testing.assert_allclose(from_ndarray, from_axis)
