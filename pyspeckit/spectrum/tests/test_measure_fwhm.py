"""
Regression tests for issue #141: measure_approximate_fwhm returned wrong
values for off-center lines because the region-growing logic mirrored the
first above-threshold pixel across the array midpoint instead of finding the
last above-threshold pixel (and the half-open slices never actually grew the
region rightward).
"""
import warnings

import numpy as np

import pyspeckit

SIGMA_TO_FWHM = 2 * np.sqrt(2 * np.log(2))


def make_gaussian_spectrum(center, sigma, npix=101, xmax=100., amp=1.,
                           noise=0., seed=None):
    xarr = np.linspace(0, xmax, npix)
    data = amp * np.exp(-(xarr - center)**2 / (2. * sigma**2))
    if noise:
        data = data + np.random.RandomState(seed).randn(data.size) * noise
    with warnings.catch_warnings():
        # ignore warning about creating an empty header
        warnings.simplefilter('ignore')
        sp = pyspeckit.Spectrum(xarr=xarr, data=data,
                                error=np.ones_like(data) * max(noise, 1e-6),
                                xarrkwargs={'unit': 'km/s'})
    return sp


def fit_and_measure(sp, center, sigma, **kwargs):
    sp.specfit(fittype='gaussian', guesses=[1., center, sigma])
    return sp.specfit.measure_approximate_fwhm(**kwargs).value


def test_fwhm_offcenter_grow_region():
    """
    Off-center narrow line with a threshold above the half-max so only one
    pixel exceeds it: the grow-by-one-pixel branch must expand around the
    *line*, not around its mirror image (issue #141).
    """
    sigma = 0.8
    center = 25.  # 25% of the axis: well off-center
    analytic = SIGMA_TO_FWHM * sigma

    sp = make_gaussian_spectrum(center, sigma)

    # interpolated measurement: quantized at 0.1 km/s, so ~5% is attainable;
    # the pre-fix bug gave ~25.9 (>1200% error)
    fwhm = fit_and_measure(sp, center, sigma, threshold=0.6,
                           interpolate_factor=10)
    assert abs(fwhm - analytic) / analytic < 0.07

    # non-interpolated measurement: quantized at the 1 km/s channel width
    fwhm_coarse = fit_and_measure(sp, center, sigma, threshold=0.6)
    assert abs(fwhm_coarse - analytic) / analytic < 0.10


def test_fwhm_offcenter_accuracy():
    """Well-resolved off-center line, no region growing: FWHM should match
    the analytic value closely with and without interpolation."""
    sigma = 5.
    center = 25.
    analytic = SIGMA_TO_FWHM * sigma

    sp = make_gaussian_spectrum(center, sigma, npix=1001)

    fwhm = fit_and_measure(sp, center, sigma, threshold=0.05,
                           interpolate_factor=10)
    assert abs(fwhm - analytic) / analytic < 0.01

    fwhm_coarse = fit_and_measure(sp, center, sigma, threshold=0.05)
    assert abs(fwhm_coarse - analytic) / analytic < 0.03


def test_fwhm_monte_carlo():
    """
    The usage pattern from issue #141: repeated noise realizations must give
    a distribution of FWHMs centered on the analytic value (before the fix,
    every realization returned a wildly wrong value).
    """
    sigma = 0.8
    center = 25.
    analytic = SIGMA_TO_FWHM * sigma

    fwhms = []
    for seed in range(20):
        sp = make_gaussian_spectrum(center, sigma, noise=0.02, seed=seed)
        fwhms.append(fit_and_measure(sp, center, sigma, threshold=0.6,
                                     interpolate_factor=10))
    fwhms = np.array(fwhms)

    assert abs(fwhms.mean() - analytic) / analytic < 0.10
    # a genuine (finite-width) distribution, not one constant wrong value
    assert fwhms.std() > 0


def test_fwhm_grow_region_at_edges():
    """Region growing at the array boundaries must not wrap or go out of
    bounds."""
    sigma = 0.8
    analytic = SIGMA_TO_FWHM * sigma

    for center in (1., 99.):
        sp = make_gaussian_spectrum(center, sigma)
        fwhm = fit_and_measure(sp, center, sigma, threshold=0.6,
                               interpolate_factor=10)
        assert abs(fwhm - analytic) / analytic < 0.10
