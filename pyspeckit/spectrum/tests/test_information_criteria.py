"""
Tests for log-likelihood and information-criterion (AIC/AICc/BIC) access on
Specfit (https://github.com/pyspeckit/pyspeckit/issues/246)
"""
import warnings

import numpy as np
import pytest
from astropy.stats import (akaike_info_criterion, akaike_info_criterion_lsq,
                           bayesian_info_criterion)

from .. import Spectrum

SIGMA = 0.1


def gaussian(x, amp, cen, wid):
    return amp * np.exp(-(x - cen)**2 / (2. * wid**2))


def make_spectrum(x, y_true, sigma=SIGMA, seed=0, with_error=True):
    rng = np.random.RandomState(seed)
    data = y_true + rng.randn(x.size) * sigma
    error = np.ones(x.size) * sigma if with_error else None
    with warnings.catch_warnings():
        # ignore warning about creating an empty header
        warnings.simplefilter('ignore')
        sp = Spectrum(xarr=x, data=data, error=error)
    return sp


class TestInformationCriteria(object):

    def setup_method(self):
        # n/k = 240/3 = 80 > 40, so astropy's akaike_info_criterion applies
        # no small-sample correction here
        self.x = np.linspace(-6, 6, 240)
        self.sp_one = make_spectrum(self.x, gaussian(self.x, 1.5, 0.0, 0.5),
                                    seed=42)
        y_two = (gaussian(self.x, 1.5, -2.5, 0.5) +
                 gaussian(self.x, 1.0, 2.5, 0.4))
        self.sp_two = make_spectrum(self.x, y_two, seed=1234)

    def test_logL_hand_computed(self):
        """logL must match the hand-computed Gaussian log-likelihood"""
        sp = self.sp_one
        sp.specfit(fittype='gaussian', guesses=[1.4, 0.1, 0.6])

        n = sp.specfit.npix_fitted
        assert n == self.x.size

        residuals = sp.data - sp.specfit.get_full_model()
        chi2 = np.sum((residuals / SIGMA)**2)
        logL_hand = -0.5 * chi2 - n * np.log(np.sqrt(2 * np.pi) * SIGMA)

        np.testing.assert_allclose(sp.specfit.logL, logL_hand, rtol=1e-10)
        # sanity check: the fit is good, so chi2 ~ n
        assert 0.7 * n < chi2 < 1.3 * n

    def test_aic_bic_hand_computed_and_astropy(self):
        """AIC/AICc/BIC must match both the textbook formulas and astropy"""
        sp = self.sp_one
        sp.specfit(fittype='gaussian', guesses=[1.4, 0.1, 0.6])

        k = sp.specfit.n_free_parameters
        n = sp.specfit.npix_fitted
        logL = sp.specfit.logL
        assert k == 3

        # hand-computed
        np.testing.assert_allclose(sp.specfit.AIC, 2 * k - 2 * logL,
                                   rtol=1e-10)
        np.testing.assert_allclose(sp.specfit.AICc,
                                   2 * k - 2 * logL
                                   + 2. * k * (k + 1) / (n - k - 1),
                                   rtol=1e-10)
        np.testing.assert_allclose(sp.specfit.BIC, k * np.log(n) - 2 * logL,
                                   rtol=1e-10)

        # astropy cross-checks; n/k = 80 > 40 so this is the uncorrected AIC
        np.testing.assert_allclose(sp.specfit.AIC,
                                   akaike_info_criterion(logL, k, n),
                                   rtol=1e-10)
        np.testing.assert_allclose(sp.specfit.BIC,
                                   bayesian_info_criterion(logL, k, n),
                                   rtol=1e-10)

    def test_aicc_matches_astropy_small_sample(self):
        """
        astropy's akaike_info_criterion automatically applies the
        small-sample correction when n/k < 40; check AICc against it
        """
        x = np.linspace(-4, 4, 60)  # n/k = 20 < 40
        sp = make_spectrum(x, gaussian(x, 1.5, 0.0, 0.5), seed=7)
        sp.specfit(fittype='gaussian', guesses=[1.4, 0.1, 0.6])

        k = sp.specfit.n_free_parameters
        n = sp.specfit.npix_fitted
        assert n == 60 and k == 3

        np.testing.assert_allclose(sp.specfit.AICc,
                                   akaike_info_criterion(sp.specfit.logL, k,
                                                         n),
                                   rtol=1e-10)

    def test_aic_matches_astropy_lsq(self):
        """
        astropy's akaike_info_criterion_lsq computes AIC from the sum of
        squared residuals, i.e., for the maximum-likelihood variance
        sigma^2 = SSR/n, and drops the constant -n/2*(1 + ln(2*pi)) from the
        log-likelihood.  If we refit using exactly that sigma as the error
        spectrum, our AIC must equal the astropy lsq AIC plus that constant.
        """
        sp = self.sp_one
        sp.specfit(fittype='gaussian', guesses=[1.4, 0.1, 0.6])

        n = sp.specfit.npix_fitted
        ssr = np.sum((sp.data - sp.specfit.get_full_model())**2)

        # uniform rescaling of the errors does not move the best-fit
        # parameters, so SSR is unchanged by this refit
        sp.error[:] = np.sqrt(ssr / n)
        sp.specfit(fittype='gaussian', guesses=list(sp.specfit.parinfo.values))

        k = sp.specfit.n_free_parameters
        assert sp.specfit.npix_fitted == n
        np.testing.assert_allclose(
            np.sum((sp.data - sp.specfit.get_full_model())**2), ssr,
            rtol=1e-6)

        np.testing.assert_allclose(
            sp.specfit.AIC,
            akaike_info_criterion_lsq(ssr, k, n) + n * (1 + np.log(2*np.pi)),
            rtol=1e-6)

    def test_model_selection_two_component_data(self):
        """on genuinely 2-component data, the 2-component fit must win"""
        sp = self.sp_two
        sp.specfit(fittype='gaussian', guesses=[1.2, -2.0, 0.6])
        aic1 = sp.specfit.AIC
        aicc1 = sp.specfit.AICc
        bic1 = sp.specfit.BIC

        sp.specfit(fittype='gaussian',
                   guesses=[1.2, -2.0, 0.6, 0.8, 2.0, 0.5])
        assert sp.specfit.n_free_parameters == 6
        aic2 = sp.specfit.AIC
        aicc2 = sp.specfit.AICc
        bic2 = sp.specfit.BIC

        assert aic2 < aic1
        assert aicc2 < aicc1
        assert bic2 < bic1

    def test_model_selection_one_component_data(self):
        """on 1-component data, the extra component must lose"""
        sp = self.sp_one
        sp.specfit(fittype='gaussian', guesses=[1.4, 0.1, 0.6])
        aic1 = sp.specfit.AIC
        bic1 = sp.specfit.BIC

        sp.specfit(fittype='gaussian',
                   guesses=[1.4, 0.1, 0.6, 0.2, 3.0, 0.5])
        aic2 = sp.specfit.AIC
        bic2 = sp.specfit.BIC

        assert aic1 < aic2
        assert bic1 < bic2

    def test_fixed_parameters_reduce_k(self):
        sp = self.sp_one
        sp.specfit(fittype='gaussian', guesses=[1.4, 0.0, 0.5],
                   fixed=[False, True, False])
        assert sp.specfit.n_free_parameters == 2
        np.testing.assert_allclose(sp.specfit.AIC,
                                   2 * 2 - 2 * sp.specfit.logL, rtol=1e-10)
        np.testing.assert_allclose(
            sp.specfit.BIC,
            2 * np.log(sp.specfit.npix_fitted) - 2 * sp.specfit.logL,
            rtol=1e-10)

    def test_tied_parameters_reduce_k(self):
        sp = self.sp_two
        sp.specfit(fittype='gaussian',
                   guesses=[1.2, -2.0, 0.45, 0.8, 2.0, 0.45],
                   tied=['', '', '', '', '', 'p[2]'])
        assert sp.specfit.n_free_parameters == 5
        np.testing.assert_allclose(sp.specfit.AIC,
                                   2 * 5 - 2 * sp.specfit.logL, rtol=1e-10)

    def test_no_error_spectrum_warns(self):
        sp = make_spectrum(self.x, gaussian(self.x, 1.5, 0.0, 0.5), seed=42,
                           with_error=False)
        sp.specfit(fittype='gaussian', guesses=[1.4, 0.1, 0.6])
        with pytest.warns(UserWarning,
                          match="only.*defined up to an additive constant"):
            logL = sp.specfit.logL
        assert np.isfinite(logL)

    def test_no_fit_raises(self):
        sp = make_spectrum(self.x, gaussian(self.x, 1.5, 0.0, 0.5), seed=42)
        with pytest.raises(AttributeError):
            sp.specfit.logL
        with pytest.raises(AttributeError):
            sp.specfit.n_free_parameters
