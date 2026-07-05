"""
Tests for the parameter-constraint examples documented in
docs/parameters.rst ("Constraining fit parameters: limits, fixed, and
tied").  These mirror the doc code blocks so that the documented behavior
(issue #321) stays correct.
"""
import numpy as np

import pyspeckit


def make_gaussian_spectrum(amplitude=5., center=50., sigma=10., stddev=0.1,
                           seed=0):
    xaxis = np.linspace(-50., 150., 100)
    np.random.seed(seed)
    synth_data = amplitude * np.exp(-(xaxis - center)**2 / (2 * sigma**2))
    noise = np.random.randn(xaxis.size) * stddev
    error = stddev * np.ones_like(synth_data)
    sp = pyspeckit.Spectrum(data=synth_data + noise, error=error, xarr=xaxis,
                            xarrkwargs={'unit': 'km/s'}, unit='K')
    return sp


def test_amplitude_limit_pegs_at_boundary():
    # true amplitude is 5, but the fit is constrained to [0, 3]
    sp = make_gaussian_spectrum()
    T, F = True, False
    sp.specfit(fittype='gaussian', guesses=[1, 45, 5],
               limited=[(T, T), (F, F), (T, F)],
               limits=[(0, 3), (0, 0), (0, 0)])
    assert sp.specfit.parinfo['AMPLITUDE0'].value == 3.0
    # the center is still recovered correctly
    np.testing.assert_allclose(sp.specfit.parinfo['SHIFT0'].value, 50,
                               atol=1)


def test_fixed_parameter_stays_fixed():
    sp = make_gaussian_spectrum()
    F = False
    sp.specfit(fittype='gaussian', guesses=[4, 45, 8],
               fixed=[F, F, True])
    assert sp.specfit.parinfo['WIDTH0'].value == 8.0
    assert sp.specfit.parinfo['WIDTH0'].fixed
    assert sp.specfit.parinfo['WIDTH0'].error == 0


def test_tied_parameters_fixed_separation():
    xaxis = np.linspace(6690., 6760., 200)
    sigma = 1.5
    synth_data = (1.5 * np.exp(-(xaxis - 6716.44)**2 / (2 * sigma**2)) +
                  1.0 * np.exp(-(xaxis - 6730.82)**2 / (2 * sigma**2)))
    np.random.seed(1)
    noise = np.random.randn(xaxis.size) * 0.05
    error = 0.05 * np.ones_like(synth_data)
    sp = pyspeckit.Spectrum(data=synth_data + noise, error=error, xarr=xaxis,
                            xarrkwargs={'unit': 'angstrom'},
                            unit='erg/s/cm2/AA')
    sp.specfit(fittype='gaussian', guesses=[1, 6716, 1, 1, 6731, 1],
               tied=['', '', '', '', 'p[1]+14.38', ''])
    separation = (sp.specfit.parinfo['SHIFT1'].value -
                  sp.specfit.parinfo['SHIFT0'].value)
    np.testing.assert_allclose(separation, 14.38)


def test_modified_parinfo_reuse():
    sp = make_gaussian_spectrum()
    sp.specfit(fittype='gaussian', guesses=[4, 45, 8])
    pi = sp.specfit.parinfo
    pi['AMPLITUDE0'].limited = (True, True)
    pi['AMPLITUDE0'].limits = (0, 3)
    pi['AMPLITUDE0'].value = 2.5
    sp.specfit(parinfo=pi)
    assert sp.specfit.parinfo['AMPLITUDE0'].value == 3.0
    assert len(sp.specfit.parinfo) == 3


def test_minpars_maxpars_equivalent_to_limits():
    sp = make_gaussian_spectrum()
    sp.specfit(fittype='gaussian', guesses=[1, 45, 5],
               limitedmin=[True, False, True],
               limitedmax=[True, False, False],
               minpars=[0, 0, 0],
               maxpars=[3, 0, 0])
    assert sp.specfit.parinfo['AMPLITUDE0'].value == 3.0
    assert sp.specfit.parinfo['AMPLITUDE0'].limits == (0, 3)
    assert tuple(sp.specfit.parinfo['AMPLITUDE0'].limited) == (True, True)
