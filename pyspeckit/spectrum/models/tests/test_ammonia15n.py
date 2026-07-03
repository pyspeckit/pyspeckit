"""
Tests for the 15NH3 model (thin wrapper re-using the NH3 machinery with
swapped spectroscopic constants).
"""
import numpy as np

import pyspeckit
from pyspeckit.spectrum.models import ammonia15n_constants
from pyspeckit.spectrum.models.ammonia import ammonia
from pyspeckit.spectrum.models.ammonia15n import ammonia15n
from pyspeckit.spectrum.units import SpectroscopicAxis


def make_xarr():
    # covers both the 15NH3 (1,1) line at 22.6249295 GHz and the (2,2)
    # line at 22.6498434 GHz
    return SpectroscopicAxis(np.linspace(22.610e9, 22.660e9, 3000),
                             unit='Hz')


def test_ammonia15n_model_spectrum():
    xarr = make_xarr()
    spec = ammonia15n(xarr, trot=12, tex=6, ntot=13.5, width=0.4)
    assert np.all(np.isfinite(spec))
    assert spec.max() > 0.1
    # the wrapper must be exactly the generic model with swapped constants
    spec2 = ammonia(xarr, trot=12, tex=6, ntot=13.5, width=0.4,
                    constants=ammonia15n_constants)
    np.testing.assert_array_equal(np.asarray(spec), np.asarray(spec2))


def test_ammonia_14n_unaffected():
    # the default (no constants kwarg) must still be 14NH3
    xarr = SpectroscopicAxis(np.linspace(23.68e9, 23.71e9, 1000), unit='Hz')
    spec = ammonia(xarr, trot=12, tex=6, ntot=14.5, width=0.4)
    assert np.all(np.isfinite(spec))
    assert spec.max() > 1  # the main (1,1) group is bright at this column


def test_ammonia15n_registry_fit():
    xarr = make_xarr()
    truth = dict(trot=12, tex=6, ntot=13.8, width=0.4, xoff_v=0.0)
    spec = ammonia15n(xarr, fortho=0.0, **truth)
    noise = np.random.RandomState(42).randn(xarr.size) * 0.003
    sp = pyspeckit.Spectrum(xarr=xarr, data=spec + noise)
    sp.specfit(fittype='ammonia15n',
               guesses=[10, 5, 13.5, 0.5, 0.05, 0.0],
               fixed=[False, False, False, False, False, True])
    fitted = dict(zip(['trot', 'tex', 'ntot', 'width', 'xoff_v', 'fortho'],
                      sp.specfit.modelpars))
    np.testing.assert_allclose(fitted['tex'], truth['tex'], rtol=0.1)
    np.testing.assert_allclose(fitted['ntot'], truth['ntot'], rtol=0.05)
    np.testing.assert_allclose(fitted['width'], truth['width'], rtol=0.1)
    np.testing.assert_allclose(fitted['xoff_v'], truth['xoff_v'], atol=0.05)
