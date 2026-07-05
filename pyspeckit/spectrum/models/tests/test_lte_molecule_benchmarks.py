"""
Benchmark tests for the LTE molecule model (issue #357).

These tests validate ``pyspeckit.spectrum.models.lte_molecule`` at three
levels, all of which run offline (no network access):

1. **Analytic limits** for the CO J=1-0 transition using published CDMS
   constants (nu = 115.2712018 GHz, A_ul = 7.203e-8 s^-1, g_u = 3,
   E_u/k = 5.53 K) and a rigid-rotor partition function
   (B_0 = 57.635968 GHz):

   - in the optically thin limit, the integrated intensity of the
     generated model inverts back to the input column density through the
     standard N_upper(integrated intensity) relation
     (Mangum & Shirley 2015 eqn 82; implemented as ``nupper_of_kkms`` /
     ``ntot_of_nupper``),
   - in the optically thick limit, the peak brightness saturates at
     J_nu(T_ex) - J_nu(T_bg),
   - the opacity is exactly linear in column density,
   - the integrated opacity profile equals the analytic
     integral-tau (i.e., the Gaussian profile is normalized), and the
     cgs and astropy-units tau implementations agree.

2. **XCLASS cross-check**: the stored outputs of the CO benchmark notebook

       https://github.com/pyspeckit/pyspeckit-tests/blob/master/xclass-benchmark/CO_Benchmark.ipynb

   (XCLASS run of 2021-09-05) recorded, for each of the lowest 12 CO
   rotational lines, the sum over 0.1-MHz channels of the XCLASS
   brightness-temperature model and of the equivalent pyspeckit model.
   Those values are stored in ``data/co_xclass_benchmark.ecsv``.  Here the
   pyspeckit model is regenerated from the CDMS catalog file
   ``data/co.cat`` (SPCAT format, CDMS tag 028503, CO v=0) and partition
   function ``data/co.qpart`` and compared against both the stored XCLASS
   sums and the stored pyspeckit sums.  The notebook found
   |XCLASS - pyspeckit| / pyspeckit <= 6.4e-5 per window; the offline
   reconstruction here reproduces the stored pyspeckit sums to <= 6.3e-5
   and the XCLASS sums to <= 1.2e-4 (the extra ~6e-5 comes from
   partition-function table rounding in co.qpart vs. the CDMS species
   table used by the notebook, and from reconstructing the XCLASS
   frequency grid).

3. **molsim cross-check** (https://github.com/bmcguir2/molsim): the same
   notebook stored the ratio of molsim's per-line opacity and Einstein-A
   values to pyspeckit's.  molsim normalizes the peak opacity by the FWHM
   linewidth where pyspeckit normalizes the Gaussian profile by
   sqrt(2 pi) sigma, so the opacities differ by exactly
   sqrt(2 pi / (8 ln 2)) = 1.0644670...; after removing that documented
   convention factor the codes agree to ~1.5e-4 in tau and ~3e-5 in A_ul.
   (A live molsim comparison is not run here: molsim is not
   pip-installable -- the PyPI name "molsim" belongs to an unrelated,
   renamed package -- so the stored notebook values are used instead.)

All reference data files live in ``data/`` next to this file; see the
ECSV header of ``co_xclass_benchmark.ecsv`` for full provenance and the
physical parameters of the benchmark model (T=79.84 K, N=2.106e17 cm^-2,
FWHM=4.21 km/s, v_cen=-1.1 km/s, T_bg=2.73 K).
"""
import os

import numpy as np
import pytest
from astropy import constants
from astropy import units as u
from astropy.table import Table

from ..lte_molecule import (Jnu, Jnu_cgs, generate_model, line_tau,
                            line_tau_cgs, ntot_of_nupper, nupper_of_kkms)

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

# np.trapz was deprecated & renamed to np.trapezoid in numpy 2.0
trapezoid = getattr(np, 'trapezoid', None) or np.trapz

CKMS = constants.c.to(u.km / u.s).value

# ---------------------------------------------------------------------------
# CO J=1-0 constants (CDMS, as quoted in issue #357)
# ---------------------------------------------------------------------------
NU0 = 115.2712018 * u.GHz
AUL = 7.203e-8 * u.Hz       # Einstein A_ul
GU = 3                      # upper-state degeneracy (2J+1)
EU_K = 5.53 * u.K           # E_u / k_B
B0 = 57.635968 * u.GHz      # rigid-rotor rotation constant
EU_ERG = (EU_K * constants.k_B).to(u.erg)


def rigid_rotor_partfunc(tem):
    """Q(T) = sum_J (2J+1) exp(-h B J(J+1) / k T) for a linear rigid rotor."""
    tem = u.Quantity(tem, u.K)
    J = np.arange(200)
    return np.sum((2 * J + 1) *
                  np.exp((-constants.h * B0 * J * (J + 1) /
                          (constants.k_B * tem)).decompose().value))


def co10_model(xarr, column, tex, sigma_kms=1.0, vcen=0.0, tbg=2.73,
               **kwargs):
    """Single-line CO 1-0 LTE model built from the constants above."""
    return generate_model(xarr, vcen, sigma_kms, tex, column,
                          freqs=u.Quantity([NU0]),
                          aij=np.log10([AUL.value]),
                          deg=np.array([GU]),
                          EU=np.array([EU_ERG.value]),
                          partfunc=rigid_rotor_partfunc,
                          tbg=tbg, **kwargs)


def co10_grid(nchan=20001, halfwidth_MHz=10.):
    """A finely-sampled frequency grid centered on the CO 1-0 line."""
    nu0 = NU0.to(u.MHz).value
    return np.linspace(nu0 - halfwidth_MHz, nu0 + halfwidth_MHz,
                       nchan) * u.MHz


# ---------------------------------------------------------------------------
# Layer 1: analytic-limit benchmarks
# ---------------------------------------------------------------------------
def test_thin_limit_column_roundtrip():
    """
    Optically thin limit: the integrated intensity of the generated model
    must invert back to the input total column density through the
    standard optically-thin column-density relation
    N_u = (8 pi nu k_B) / (A_ul c^2 h) * int(T_R dnu)
    (Mangum & Shirley 2015 eqn 82, implemented in nupper_of_kkms) followed
    by the Boltzmann/partition-function correction (eqn 31,
    ntot_of_nupper).

    tbg=0.1 K makes J_nu(T_bg) ~ 1e-23 K, so the negligible-background
    assumption of the relation holds to machine precision.  The measured
    round-trip error is ~5e-6 (dominated by the finite peak opacity,
    tau_0 ~ 1.3e-5); assert the <1% requirement with lots of margin.
    """
    column = 1e12  # cm^-2; gives peak tau ~ 1.3e-5 at 25 K
    tex = 25.0
    xarr = co10_grid()
    model = co10_model(xarr, column, tex, tbg=0.1)

    # integrated intensity in K km/s: int T dnu * (c / nu0)
    int_T_dnu = trapezoid(model, xarr.to(u.Hz).value)  # K Hz
    kkms = int_T_dnu * CKMS / NU0.to(u.Hz).value * u.K * u.km / u.s

    nupper = nupper_of_kkms(kkms, NU0, AUL)
    ntot = ntot_of_nupper(nupper, EU_ERG, tex * u.K,
                          Q_rot=rigid_rotor_partfunc(tex), degeneracy=GU)

    np.testing.assert_allclose(ntot.to(u.cm ** -2).value, column, rtol=0.01)
    # the actual agreement is much better than the 1% requirement
    np.testing.assert_allclose(ntot.to(u.cm ** -2).value, column, rtol=1e-4)


def test_thin_limit_integrated_intensity_analytic():
    """
    Optically thin limit with a background: the integrated model intensity
    must equal (J_nu(T_ex) - J_nu(T_bg)) * int(tau dnu) with
    int(tau dnu) = c^2/(8 pi nu^2) A_ul N_u (exp(h nu/k T_ex) - 1)
    (Mangum & Shirley 2015 eqns 24 & 29).  Measured agreement ~5e-6;
    assert <0.1%.
    """
    column = 1e12
    tex = 25.0
    tbg = 2.73
    xarr = co10_grid()
    model = co10_model(xarr, column, tex, tbg=tbg)

    int_T_dnu = trapezoid(model, xarr.to(u.Hz).value) * u.K * u.Hz

    taudnu = line_tau(tex=tex * u.K, total_column=column * u.cm ** -2,
                      partition_function=rigid_rotor_partfunc(tex),
                      degeneracy=GU, frequency=NU0, energy_upper=EU_ERG,
                      einstein_A=AUL)
    analytic = ((Jnu(NU0, tex * u.K) - Jnu(NU0, tbg * u.K)) *
                taudnu).to(u.K * u.Hz)

    np.testing.assert_allclose(int_T_dnu.value, analytic.value, rtol=1e-3)


def test_thick_limit_peak_brightness():
    """
    High optical depth: the peak brightness must saturate at
    J_nu(T_ex) - J_nu(T_bg).  With N = 1e19 cm^-2 the peak tau is ~1e5,
    so 1 - exp(-tau) = 1 to machine precision; assert <0.1% (measured
    agreement is ~1e-9 when J_nu is evaluated at the peak channel's
    frequency).
    """
    column = 1e19
    tex = 40.0
    tbg = 2.73
    xarr = co10_grid()
    model = co10_model(xarr, column, tex, tbg=tbg)

    peak_chan = np.argmax(model)
    nu_peak = xarr[peak_chan].to(u.Hz)
    expected = (Jnu(nu_peak, tex * u.K) - Jnu(nu_peak, tbg * u.K)).to(u.K)

    np.testing.assert_allclose(model[peak_chan], expected.value, rtol=1e-3)


def test_tau_linear_in_column():
    """The stick (integrated) opacity must scale exactly linearly with N."""
    tex = 25.0
    xarr = co10_grid(nchan=101)
    tau1 = co10_model(xarr, 1e12, tex, get_tau_sticks=True)[0]
    tau10 = co10_model(xarr, 1e13, tex, get_tau_sticks=True)[0]
    tau1000 = co10_model(xarr, 1e15, tex, get_tau_sticks=True)[0]

    np.testing.assert_allclose(tau10 / tau1, 10.0, rtol=1e-10)
    np.testing.assert_allclose(tau1000 / tau1, 1000.0, rtol=1e-10)


def test_integrated_tau_consistency():
    """
    The opacity profile returned by get_tau must integrate (over
    frequency) to the stick value returned by get_tau_sticks, i.e. the
    Gaussian profile must be normalized, and both the cgs and
    astropy-units implementations of tau-dnu must agree with each other
    and with a from-scratch evaluation of MS15 eqn 29.
    """
    column = 1e14
    tex = 25.0
    xarr = co10_grid()

    tau_profile = co10_model(xarr, column, tex, get_tau=True)
    tau_stick = co10_model(xarr, column, tex, get_tau_sticks=True)[0]

    int_tau_dnu = trapezoid(tau_profile, xarr.to(u.Hz).value)
    np.testing.assert_allclose(int_tau_dnu, tau_stick, rtol=1e-6)

    # cross-check the two implementations of the same equation
    Q = rigid_rotor_partfunc(tex)
    taudnu_units = line_tau(tex=tex * u.K, total_column=column * u.cm ** -2,
                            partition_function=Q, degeneracy=GU,
                            frequency=NU0, energy_upper=EU_ERG,
                            einstein_A=AUL).to(u.Hz)
    taudnu_cgs = line_tau_cgs(tex=tex, total_column=column,
                              partition_function=Q, degeneracy=GU,
                              frequency=NU0.to(u.Hz).value,
                              energy_upper=EU_ERG.value,
                              einstein_A=AUL.value)
    np.testing.assert_allclose(taudnu_units.value, taudnu_cgs, rtol=1e-8)
    np.testing.assert_allclose(tau_stick, taudnu_cgs, rtol=1e-8)

    # from-scratch: int tau dnu = c^2/(8 pi nu^2) A N_u (e^(h nu/kT) - 1)
    nu = NU0.to(u.Hz).value
    N_u = (column * GU / Q *
           np.exp(-(EU_ERG / (constants.k_B * tex * u.K)).decompose().value))
    c_cgs = constants.c.cgs.value
    hoverk = (constants.h / constants.k_B).cgs.value
    byhand = (c_cgs ** 2 / (8 * np.pi * nu ** 2) * AUL.value * N_u *
              (np.exp(hoverk * nu / tex) - 1))
    np.testing.assert_allclose(taudnu_cgs, byhand, rtol=1e-10)


def test_jnu_implementations_agree():
    """Jnu (astropy-units) and Jnu_cgs must agree."""
    for T in (2.73, 25., 300.):
        np.testing.assert_allclose(
            Jnu(NU0, T * u.K).to(u.K).value,
            Jnu_cgs(NU0.to(u.Hz).value, T), rtol=1e-10)


# ---------------------------------------------------------------------------
# Layer 2 & 3 fixtures: offline CO catalog + stored XCLASS/molsim references
# ---------------------------------------------------------------------------
def read_spcat(catfile):
    """
    Parse an SPCAT-format .cat file (fixed-width; Pickett et al. 1998)
    and derive (freqs, log10(Aij), gup, E_upper[erg], partfunc-input
    columns) exactly the way
    ``lte_molecule.get_molecular_parameters`` does from the equivalent
    astroquery columns.
    """
    freqs, lgint, elo_icm, gup = [], [], [], []
    with open(catfile, 'r') as fh:
        for line in fh:
            if not line.strip():
                continue
            freqs.append(float(line[0:13]))    # FREQ [MHz]
            lgint.append(float(line[21:29]))   # LGINT [log10(nm^2 MHz)]
            elo_icm.append(float(line[31:41]))  # ELO [cm^-1]
            gup.append(int(line[41:44]))       # GUP
    return (np.array(freqs) * u.MHz, np.array(lgint), np.array(elo_icm),
            np.array(gup, dtype=float))


def cdms_aij_eu(freqs, lgint, elo_icm, gup, partfunc):
    """
    Convert CDMS/JPL LGINT to log10(Einstein A) and upper-state energy,
    replicating get_molecular_parameters (which follows the CDMS
    documentation and B. McGuire's simulate_lte).
    """
    freq_MHz = freqs.to(u.MHz).value
    eupper_icm = elo_icm + freqs.to(u.cm ** -1, u.spectral()).value
    EU = (eupper_icm / u.cm).to(u.erg, u.spectral()).value
    CT = 300
    sijmu = ((np.exp(-(elo_icm / 0.695) / CT) -
              np.exp(-(eupper_icm / 0.695) / CT)) ** (-1) *
             ((10 ** lgint) / freq_MHz) * 24025.120666 * partfunc(CT))
    aij = np.log10(1.16395e-20 * freq_MHz ** 3 * sijmu / gup)
    return aij, EU


@pytest.fixture(scope='module')
def co_benchmark():
    """
    Load the CO catalog + reference table and regenerate the pyspeckit
    model of the benchmark notebook on each of the 12 spectral windows.
    """
    reftab = Table.read(os.path.join(DATA_DIR, 'co_xclass_benchmark.ecsv'))
    pars = reftab.meta['parameters']

    freqs, lgint, elo_icm, gup = read_spcat(os.path.join(DATA_DIR, 'co.cat'))

    qtab = np.loadtxt(os.path.join(DATA_DIR, 'co.qpart'), comments='#')
    qtems, qvals = qtab[:, 0], qtab[:, 1]

    def partfunc(tem):
        # linear interpolation in Q vs T, matching the partfunc returned
        # by get_molecular_parameters
        return np.interp(u.Quantity(tem, u.K).value, qtems, qvals)

    aij, EU = cdms_aij_eu(freqs, lgint, elo_icm, gup, partfunc)

    tkin = pars['tkin_K']
    ntot = pars['Ntot_cm-2']
    sigma = pars['vwidth_fwhm_kms'] / np.sqrt(8 * np.log(2))
    vcen = pars['vcen_kms']
    tbg = pars['tbg_K']
    step = pars['freqstep_MHz']

    sums = []
    peak_taus = []
    for row in reftab:
        grid = np.arange(row['fmin_MHz'], row['fmax_MHz'] + step / 2,
                         step) * u.MHz
        model = generate_model(grid, vcen, sigma, tkin, ntot, freqs=freqs,
                               aij=aij, deg=gup, EU=EU, partfunc=partfunc,
                               tbg=tbg)
        tau = generate_model(grid, vcen, sigma, tkin, ntot, freqs=freqs,
                             aij=aij, deg=gup, EU=EU, partfunc=partfunc,
                             tbg=tbg, get_tau=True)
        sums.append(model.sum())
        peak_taus.append(tau.max())

    return {'reftab': reftab, 'freqs': freqs, 'gup': gup, 'aij': aij,
            'EU': EU, 'partfunc': partfunc,
            'sums': np.array(sums), 'peak_taus': np.array(peak_taus)}


# ---------------------------------------------------------------------------
# Layer 2: XCLASS cross-check
# ---------------------------------------------------------------------------
def test_catalog_constants_match_co10(co_benchmark):
    """The shipped co.cat must reproduce the published CO 1-0 constants."""
    freqs = co_benchmark['freqs']
    assert freqs[0].to(u.MHz).value == 115271.2018
    assert co_benchmark['gup'][0] == GU
    np.testing.assert_allclose(10 ** co_benchmark['aij'][0], AUL.value,
                               rtol=1e-3)
    eu_K = (co_benchmark['EU'][0] * u.erg / constants.k_B).to(u.K).value
    np.testing.assert_allclose(eu_K, EU_K.value, rtol=1e-3)


def test_catalog_frequencies_match_reference(co_benchmark):
    """The 12 benchmark lines must be the first 12 entries of co.cat."""
    np.testing.assert_allclose(
        co_benchmark['freqs'][:12].to(u.MHz).value,
        co_benchmark['reftab']['restfreq_MHz'], rtol=0, atol=1e-4)


def test_reproduces_stored_pyspeckit_sums(co_benchmark):
    """
    Regression: the offline reconstruction (co.cat + co.qpart) must
    reproduce the pyspeckit-CDMS window sums stored by the benchmark
    notebook.  Measured agreement: <= 6.3e-5 (limited by co.qpart
    partition-function rounding relative to the CDMS species table the
    notebook queried).
    """
    np.testing.assert_allclose(co_benchmark['sums'],
                               co_benchmark['reftab']['pyspeckit_cdms_sum_K'],
                               rtol=1.5e-4)


def test_agrees_with_xclass(co_benchmark):
    """
    Cross-code benchmark: the pyspeckit LTE model must match the stored
    XCLASS brightness-temperature window sums.  The notebook measured
    |XCLASS - pyspeckit|/pyspeckit <= 6.4e-5 on the native XCLASS grid;
    the offline reconstruction here achieves <= 1.2e-4.
    """
    np.testing.assert_allclose(co_benchmark['sums'],
                               co_benchmark['reftab']['xclass_sum_K'],
                               rtol=3e-4)


def test_reproduces_stored_peak_taus(co_benchmark):
    """
    The recomputed peak channel opacities must match the values stored by
    the notebook.  The stored values were sampled on the full-band molsim
    grid while these are sampled on per-window grids, so channel-center
    phase differences limit the agreement to ~1e-3 (measured <= 9e-4).
    The notebook also recorded that the XCLASS per-line opacities match
    pyspeckit's to within [-2.9e-4, +1.6e-4] (ratio histogram in the
    notebook), i.e. the tau agreement with XCLASS is at the few-1e-4
    level.
    """
    np.testing.assert_allclose(co_benchmark['peak_taus'],
                               co_benchmark['reftab']['pyspeckit_peak_tau'],
                               rtol=2.5e-3)


# ---------------------------------------------------------------------------
# Layer 3: molsim cross-check (stored values; see module docstring for why
# a live molsim run is not attempted)
# ---------------------------------------------------------------------------
MOLSIM_TAU_CONVENTION = np.sqrt(2 * np.pi / (8 * np.log(2)))


def test_molsim_tau_convention_factor(co_benchmark):
    """
    molsim's stored opacities differ from pyspeckit's by exactly the
    profile-normalization convention factor sqrt(2 pi/(8 ln 2)) (molsim
    divides by the FWHM, pyspeckit by sqrt(2 pi) sigma).  After removing
    it, the stored ratios agree to ~1.5e-4.
    """
    ratio = co_benchmark['reftab']['molsim_tau_ratio']
    np.testing.assert_allclose(ratio, MOLSIM_TAU_CONVENTION, rtol=3e-4)

    # tie the stored molsim opacities to the live pyspeckit computation:
    molsim_tau = (co_benchmark['reftab']['molsim_tau_ratio'] *
                  co_benchmark['reftab']['pyspeckit_peak_tau'])
    np.testing.assert_allclose(co_benchmark['peak_taus'] *
                               MOLSIM_TAU_CONVENTION,
                               molsim_tau, rtol=2.5e-3)


def test_molsim_einstein_A_agreement(co_benchmark):
    """
    molsim's Einstein A values (read directly from the SPCAT catalog) must
    agree with the ones pyspeckit derives from the CDMS LGINT column.
    The stored ratios agree to <= 3.2e-5; assert <= 1e-4.
    """
    ratio = co_benchmark['reftab']['molsim_aij_ratio']
    np.testing.assert_allclose(ratio, 1.0, rtol=1e-4)
