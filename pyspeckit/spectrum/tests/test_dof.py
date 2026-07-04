"""
Regression tests for issue #258: the fitter could report negative degrees of
freedom.

The root cause: ``npix_fitted`` was computed as ``includemask.sum() -
mask.sum()``, where ``mask`` (non-finite data pixels) was summed over the
*whole* spectrum.  A spectrum with more blanked (NaN) pixels outside the
fitted region than pixels inside it went negative, silently producing a
negative dof and negative reduced chi^2.

Separately, a fit whose included region genuinely contains fewer pixels than
the model has free parameters is under-determined; dof <= 0 now emits a
warning (but the arithmetic is not clamped).
"""
import warnings

import numpy as np
import pytest

import pyspeckit


def make_two_gaussian_spectrum(nan_edge=0):
    xarr = pyspeckit.units.SpectroscopicAxis(np.linspace(-25, 25, 100),
                                             unit='km/s')
    data = (5 * np.exp(-xarr.value**2 / 2.) +
            3 * np.exp(-(xarr.value - 5)**2 / 2.))
    if nan_edge:
        data[:nan_edge] = np.nan
    with warnings.catch_warnings():
        # ignore warning about creating an empty header
        warnings.simplefilter('ignore')
        return pyspeckit.Spectrum(data=data, xarr=xarr,
                                  error=np.ones(data.size) * 0.1)


GUESSES_2G = [5, 0, 1, 3, 5, 1]  # 2 gaussians = 6 parameters


def test_dof_ignores_masked_pixels_outside_selection():
    # Regression test for issue #258: 20 NaN pixels at the band edge, well
    # outside the fitted region, used to be subtracted from npix_fitted,
    # yielding npix_fitted = 14 - 20 = -6 and dof = -12.
    sp = make_two_gaussian_spectrum(nan_edge=20)
    with warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter('always')
        sp.specfit(fittype='gaussian', guesses=GUESSES_2G,
                   xmin=-1.5, xmax=5.5)
        dof = sp.specfit.dof
    assert sp.specfit.includemask.sum() == 14
    assert sp.specfit.mask.sum() == 20
    assert sp.specfit.npix_fitted == 14
    assert dof == 14 - 6
    assert not any('under-determined' in str(w.message) for w in wlist)


def test_dof_counts_masked_pixels_inside_selection():
    # NaN pixels *inside* the included region do reduce npix_fitted
    sp = make_two_gaussian_spectrum()
    sp.data[50] = np.nan  # x = +0.25 km/s, inside the fit region
    sp.specfit(fittype='gaussian', guesses=GUESSES_2G, xmin=-5, xmax=10)
    included = sp.specfit.includemask.sum()
    assert sp.specfit.npix_fitted == included - 1
    with warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter('always')
        assert sp.specfit.dof == included - 1 - 6
        assert not any('under-determined' in str(w.message) for w in wlist)


def test_underdetermined_fit_warns():
    # A wide window with most of it excluded leaves only 4 included pixels
    # for a 6-parameter model: dof = 4 - 6 = -2.  The dof arithmetic must be
    # preserved (no silent clamping), but a warning must be emitted.
    sp = make_two_gaussian_spectrum()
    with pytest.warns(UserWarning, match='under-determined'):
        sp.specfit(fittype='gaussian', guesses=GUESSES_2G,
                   xmin=-2., xmax=12., exclude=[-1, 11])
    assert sp.specfit.includemask.sum() == 4
    assert sp.specfit.npix_fitted == 4
    assert sp.specfit.n_free_parameters == 6
    with pytest.warns(UserWarning, match='4 pixels.*6 free parameters'):
        assert sp.specfit.dof == 4 - 6


def test_dof_normal_fit_no_warning():
    sp = make_two_gaussian_spectrum()
    with warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter('always')
        sp.specfit(fittype='gaussian', guesses=GUESSES_2G)
        n = sp.specfit.npix_fitted
        assert n == sp.specfit.includemask.sum() == 100
        assert sp.specfit.n_free_parameters == 6
        assert sp.specfit.dof == n - 6
        assert not any('under-determined' in str(w.message) for w in wlist)


def test_dof_respects_fixed_and_tied():
    sp = make_two_gaussian_spectrum()
    # fix the width of the first component; tie the second amplitude to the
    # first => 6 parameters, only 4 free
    sp.specfit(fittype='gaussian', guesses=GUESSES_2G,
               fixed=[False, False, True, False, False, False],
               tied=['', '', '', '0.6*p[0]', '', ''])
    assert sp.specfit.n_free_parameters == 4
    assert sp.specfit.dof == sp.specfit.npix_fitted - 4


def test_optimal_chi2_dof_counts_all_free_parameters():
    # optimal_chi2 used self.fitter.npars (per-component: 3 for a gaussian),
    # ignoring npeaks, fixed, and tied parameters
    sp = make_two_gaussian_spectrum()
    sp.specfit(fittype='gaussian', guesses=GUESSES_2G)
    modelmask = sp.specfit._compare_to_threshold(threshold='error')
    chi2 = sp.specfit.optimal_chi2(reduced=False, threshold='error')
    reduced = sp.specfit.optimal_chi2(reduced=True, threshold='error')
    expected_dof = modelmask.sum() - 6
    assert expected_dof > 0
    np.testing.assert_allclose(reduced, chi2 / expected_dof)
