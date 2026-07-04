"""
Regression tests for a batch of long-standing crash bugs in the 1-D spectrum
core (fixed together; each test names the file it guards).
"""
import warnings

import numpy as np
import pytest
from astropy.io import fits

import pyspeckit
from pyspeckit.spectrum.classes import Spectra, Spectrum
from pyspeckit.spectrum import headers


def make_gaussian_spectrum():
    xarr = pyspeckit.units.SpectroscopicAxis(np.linspace(-50, 50, 101),
                                             unit='km/s',
                                             refX=23.7e9, refX_unit='Hz')
    data = (np.exp(-np.linspace(-50, 50, 101)**2 / 50.) +
            np.random.RandomState(42).randn(101) * 0.01)
    return pyspeckit.Spectrum(xarr=xarr, data=data)


def test_multifit_keyword_deprecation():
    # fitters.py: log.warning(msg, DeprecationWarning) crashed astropy's logger
    sp = make_gaussian_spectrum()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        sp.specfit(fittype='gaussian', guesses=[1, 0, 5], multifit=True)
    assert any(issubclass(wi.category, DeprecationWarning) for wi in w)


def test_npeaks_is_integer():
    # fitters.py: py2 true-division made npeaks a float
    sp = make_gaussian_spectrum()
    sp.specfit(fittype='gaussian', guesses=[1, 0, 5])
    assert sp.specfit.npeaks == 1
    assert isinstance(sp.specfit.npeaks, (int, np.integer))


def test_get_model_xlimits_unit():
    # fitters.py: get_model_xlimits used the deprecated (default-None)
    # ``units`` variable instead of ``unit`` for the conversion
    sp = make_gaussian_spectrum()
    sp.specfit(fittype='gaussian', guesses=[1, 0, 5])
    limits = sp.specfit.get_model_xlimits(unit='km/s')
    assert limits.unit.is_equivalent('km/s')


def test_interpnans():
    # interpolation.py: interpnans used py2-era boolean subtraction and
    # could reference OK before assignment
    sp = make_gaussian_spectrum()
    sp.data[10] = np.nan
    sp.data[50] = np.inf
    sp.interpnans()
    assert np.all(np.isfinite(sp.data))


def test_interpnans_masked():
    sp = make_gaussian_spectrum()
    sp.data = np.ma.masked_array(sp.data, mask=np.zeros(sp.data.size,
                                                        dtype=bool))
    sp.data.mask[3] = True
    sp.interpnans()
    assert np.all(np.isfinite(sp.data))


def test_spline_baseline():
    # baseline.py: py2 true-division produced a float slice index in _spline
    sp = make_gaussian_spectrum()
    sp.baseline(spline=True, order=3, spline_sampling=10)
    assert sp.baseline.subtracted


def test_spectra_add_spectrum():
    # classes.py: Spectra.__add__ had an or-precedence bug making the
    # Spectrum branch unreachable (and returned None)
    sp1 = make_gaussian_spectrum()
    sp2 = make_gaussian_spectrum()
    sp3 = make_gaussian_spectrum()
    spectra = Spectra([sp1, sp2])
    result = spectra + sp3
    assert result is not None
    assert len(result.speclist) == 3
    assert isinstance(result.speclist[-1], Spectrum)


def test_as_unit_center_frequency_unit():
    # units.py: inverted branch discarded the user-supplied
    # center_frequency_unit
    sp = make_gaussian_spectrum()
    xnew = sp.xarr.as_unit('GHz', center_frequency=23.7,
                           center_frequency_unit='GHz',
                           velocity_convention='radio')
    assert xnew.unit.is_equivalent('GHz')
    np.testing.assert_allclose(xnew.value.mean(), 23.7, rtol=1e-3)


@pytest.mark.parametrize('convention', ['radio', 'optical', 'relativistic'])
def test_frequency_to_velocity_conventions(convention):
    # units.py: the relativistic-convention formula had a wrong denominator,
    # (f0**2 + f)**2 instead of (f0**2 + f**2); cross-check every convention
    # against astropy's doppler equivalencies
    from astropy import units as u
    from pyspeckit.spectrum.units import (frequency_to_velocity,
                                          velocity_conventions)
    restfreq = 23.7 * u.GHz
    freqs = np.linspace(23.5, 23.9, 11)
    expected = (freqs * u.GHz).to(
        u.km / u.s, equivalencies=velocity_conventions[convention](restfreq))
    result = frequency_to_velocity(freqs, 'GHz',
                                   center_frequency=restfreq.value,
                                   center_frequency_units='GHz',
                                   velocity_units='km/s',
                                   convention=convention)
    np.testing.assert_allclose(np.asarray(result), expected.value,
                               rtol=1e-10, atol=1e-6)


def test_polynomial_continuum_fit():
    # Regression test for issue #334: polynomial_continuum.polymodel called
    # numpy.polyval directly on a SpectroscopicAxis, which raises
    # "no implementation found for 'numpy.polyval'" because polyval is not
    # supported through __array_function__ for Quantity subclasses.
    # This is the documented workflow from
    # docs/example_continuum_fromscratch.py ("Fitting a continuum model as a
    # model").
    xaxis = np.linspace(-50, 150, 100)
    sigma = 10.
    center = 50.
    baseline = np.poly1d([0.1, 0.25])(xaxis)
    synth_data = np.exp(-(xaxis - center)**2 / (sigma**2 * 2.)) + baseline
    stddev = 0.1
    noise = np.random.RandomState(42).randn(xaxis.size) * stddev
    error = stddev * np.ones_like(synth_data)
    data = noise + synth_data

    sp = pyspeckit.Spectrum(data=data, error=error, xarr=xaxis,
                            xarrkwargs={'unit': 'km/s'},
                            unit='erg/s/cm^2/AA')
    sp.specfit.Registry.add_fitter(
        'polycontinuum',
        pyspeckit.models.polynomial_continuum.poly_fitter(),
        2)
    sp.specfit(fittype='polycontinuum', guesses=(0, 0), exclude=[30, 70])

    # recovered coefficients should be close to the input [0.1, 0.25]
    np.testing.assert_allclose(sp.specfit.parinfo.values, [0.1, 0.25],
                               atol=0.05)

    # continuum subtraction (the second half of the documented workflow).
    # The copy inherits the exclude=[30, 70] region selection; reset=True
    # restores the full range before fitting the line that lives inside the
    # previously-excluded window (the docs example achieves the same via the
    # plotter's displayed limits).
    sp_contsub = sp.copy()
    sp_contsub.data -= sp.specfit.get_full_model()
    sp_contsub.specfit(fittype='gaussian', guesses=[1, 50, 10], reset=True)
    np.testing.assert_allclose(sp_contsub.specfit.parinfo.values,
                               [1., 50., 10.], rtol=0.2)


def test_headers_intersection_conflict():
    # headers.py: NameError (Header2 vs header2) in the if_conflict=2 branch
    h1 = fits.Header()
    h2 = fits.Header()
    h1['CRVAL1'] = 1.0
    h2['CRVAL1'] = 2.0
    newheader = headers.intersection(h1, h2, if_conflict=2)
    assert newheader['CRVAL1'] == 2.0
