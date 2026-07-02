"""
Tests for the RADEX-grid ammonia model (ammonia_grids / ammonia_radex /
ammonia_model_radex), using a synthetic analytic model grid.
"""
import itertools

import numpy as np
import pytest
from astropy import units as u
from astropy.table import Table

from ...classes import Spectrum, units
from ..ammonia import ammonia_radex, ammonia_model_radex
from ..ammonia_grids import ammonia_grids, parbounds, TCMB

LINE_IDS = ['11', '22', '33', '44', '55', '66']
PARA_IDS = ['11', '22', '44', '55']
ORTHO_IDS = ['33', '66']


def tau_func(ii, logn, logT, logN, logs):
    """Smooth analytic opacity as a function of the (log) grid parameters"""
    return (0.1 * (logN - 12) + 0.01 * logn + 0.002 * ii
            - 0.02 * logs + 0.005 * logT)


def tex_func(ii, logn, logT, logN, logs):
    """Smooth analytic Tex as a function of the (log) grid parameters"""
    return (5.0 + (logN - 13) + 0.1 * logn + 0.05 * ii
            + 0.3 * logT + 0.1 * logs)


@pytest.fixture(scope='module')
def gridfile(tmp_path_factory):
    """
    Build a synthetic RADEX grid file: the first 12 columns are the model
    outputs (tau_11...tau_66, Tex_11...Tex_66); the last four columns are
    the model parameters (nH2, Temperature, Column, sigmav), as expected by
    ammonia_grids.
    """
    nh2 = np.logspace(3, 6, 4)
    temperature = np.array([5., 10., 20., 40.])
    column = np.logspace(13, 15, 5)
    sigmav = np.logspace(-1, 0.5, 4)

    prod = np.array(list(itertools.product(nh2, temperature, column, sigmav)))
    nn, TT, NN, ss = prod.T
    logn, logT, logN, logs = np.log10(nn), np.log10(TT), np.log10(NN), np.log10(ss)

    names = (['tau_' + lid for lid in LINE_IDS] +
             ['Tex_' + lid for lid in LINE_IDS] +
             ['nH2', 'Temperature', 'Column', 'sigmav'])
    columns = ([tau_func(ii + 1, logn, logT, logN, logs)
                for ii in range(len(LINE_IDS))] +
               [tex_func(ii + 1, logn, logT, logN, logs)
                for ii in range(len(LINE_IDS))] +
               [nn, TT, NN, ss])

    tbl = Table(columns, names=names)
    fn = str(tmp_path_factory.mktemp('radexgrid') / 'grid.fits')
    tbl.write(fn)
    return fn


@pytest.fixture(scope='module')
def sampler(gridfile):
    samp = ammonia_grids(gridfile=gridfile)
    assert samp is not None
    return samp


def test_sampler_exact_at_gridpoint(sampler):
    """
    Order-1 interpolation must be exact at grid nodes.
    """
    logdens, tkin, ntot, sigma = 4.0, 10.0, 14.0, 1.0
    out = sampler(logdens=logdens, tkin=tkin, ntot=ntot, fortho=0.0,
                  sigma=sigma, scaling_method='none')

    logn, logT, logN, logs = logdens, np.log10(tkin), ntot, np.log10(sigma)
    for jj, lid in enumerate(PARA_IDS):
        ii = LINE_IDS.index(lid) + 1
        np.testing.assert_allclose(out['tau_' + lid],
                                   tau_func(ii, logn, logT, logN, logs),
                                   rtol=1e-6)
        np.testing.assert_allclose(out['Tex_' + lid],
                                   tex_func(ii, logn, logT, logN, logs),
                                   rtol=1e-6)
    # fortho = 0: ortho lines have no emission
    for lid in ORTHO_IDS:
        assert out['tau_' + lid] == 0.0
        assert out['Tex_' + lid] == TCMB


def test_fortho_invariance_para(sampler):
    """
    The para-line outputs must depend only on the para column
    N_para = 10**ntot * (1 - fortho), independent of fortho.
    (regression test for the fortho bug reported in PR #332)
    """
    N_para = 10**13.8
    reference = None
    for fortho in (0.1, 0.5, 0.9):
        ntot = np.log10(N_para / (1 - fortho))
        out = sampler(logdens=4.2, tkin=13., ntot=ntot, fortho=fortho,
                      sigma=0.7, scaling_method='none')
        vals = np.array([out[key] for key in
                         ('tau_11', 'tau_22', 'Tex_11', 'Tex_22')])
        assert np.all(np.isfinite(vals))
        if reference is None:
            reference = vals
        else:
            np.testing.assert_allclose(vals, reference, rtol=1e-6)


def test_fortho_invariance_ortho(sampler):
    """
    Same as above for the ortho lines: they must depend only on
    N_ortho = 10**ntot * fortho.
    """
    N_ortho = 10**13.8
    reference = None
    for fortho in (0.1, 0.5, 0.9):
        ntot = np.log10(N_ortho / fortho)
        out = sampler(logdens=4.2, tkin=13., ntot=ntot, fortho=fortho,
                      sigma=0.7, scaling_method='none')
        vals = np.array([out[key] for key in
                         ('tau_33', 'tau_66', 'Tex_33', 'Tex_66')])
        assert np.all(np.isfinite(vals))
        if reference is None:
            reference = vals
        else:
            np.testing.assert_allclose(vals, reference, rtol=1e-6)


def test_fortho_edge_robustness(sampler):
    """
    Very small (but nonzero) fortho puts the ortho column far below the
    grid; the sampler must return finite values with tau=0 for the ortho
    lines rather than NaN.  fortho=0 and fortho=1 must also return
    all-finite dictionaries.
    """
    out = sampler(logdens=4.0, tkin=10., ntot=14., fortho=1e-8,
                  sigma=1.0, scaling_method='none')
    assert all(np.isfinite(v) for v in out.values())
    for lid in ORTHO_IDS:
        assert out['tau_' + lid] == 0.0
        assert out['Tex_' + lid] == TCMB
    # para lines unaffected by the tiny ortho fraction
    for lid in PARA_IDS:
        assert out['tau_' + lid] > 0

    for fortho in (0.0, 1.0):
        out = sampler(logdens=4.0, tkin=10., ntot=14., fortho=fortho,
                      sigma=1.0, scaling_method='none')
        assert all(np.isfinite(v) for v in out.values())
        # Tex of the skipped species must be finite and nonzero to avoid
        # division-by-zero in the radiative transfer
        skipped = ORTHO_IDS if fortho == 0 else PARA_IDS
        for lid in skipped:
            assert out['tau_' + lid] == 0.0
            assert out['Tex_' + lid] == TCMB


def test_parbounds(gridfile):
    bounds = parbounds(gridfile)
    np.testing.assert_allclose(bounds['logdens'], [3, 6])
    np.testing.assert_allclose(bounds['tkin'], [5, 40])
    np.testing.assert_allclose(bounds['ntot'], [13, 15])
    np.testing.assert_allclose(bounds['sigmav'], [0.1, 10**0.5])


def test_radex_model_fit(gridfile, sampler):
    """
    End-to-end: synthesize a spectrum from the grid model, add noise, fit
    it with ammonia_model_radex, and check that the fitted model reproduces
    the noiseless input.  (Parameter recovery can be degenerate, so we
    assert on the model spectrum, not the parameters.)
    """
    np.random.seed(42)

    # cover the NH3 (1,1) and (2,2) lines at 23.6944955 & 23.7226336 GHz
    xarr = np.linspace(23.68, 23.73, 2000) * u.GHz

    truepars = dict(tkin=15., ntot=14.2, logdens=4.5, width=0.6,
                    xoff_v=0.0, fortho=0.0)
    noiseless = ammonia_radex(units.SpectroscopicAxis(xarr),
                              interpolator=sampler, **truepars)
    assert noiseless.max() > 0.5  # sanity check: there are real lines here

    stddev = 0.01
    noise = np.random.randn(xarr.size) * stddev

    sp = Spectrum(xarr=xarr, data=noiseless + noise,
                  error=np.ones(xarr.size) * stddev)

    mod = ammonia_model_radex(gridfile=gridfile)
    sp.Registry.add_fitter('ammonia_radex', mod, 6)

    sp.specfit(fittype='ammonia_radex',
               guesses=[12., 14.0, 4.0, 0.4, 0.05, 0.0],
               fixed=[False, False, False, False, False, True])

    # the fitted model must reproduce the noiseless input spectrum
    np.testing.assert_allclose(sp.specfit.model, noiseless,
                               atol=5 * stddev)
