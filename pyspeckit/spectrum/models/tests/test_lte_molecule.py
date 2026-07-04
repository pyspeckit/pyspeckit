"""
Tests for lte_molecule.get_molecular_parameters.

The offline tests monkeypatch the astroquery query tools with small fake
species/catalog tables so that the name/tag resolution and table-handling
logic is exercised without network access (regression tests for issues
#421 and #365).

The tests marked ``remote_data`` query the live JPL/CDMS services and are
only run when pytest is invoked with ``--remote-data``.
"""
import numpy as np
import pytest

from astropy import units as u
from astropy.table import Table

astroquery = pytest.importorskip("astroquery")

from ..lte_molecule import get_molecular_parameters


JPL_TEMPERATURES = [300, 225, 150, 75, 37.5, 18.5, 9.375]


def fake_jpl_species_table():
    tab = Table()
    tab['TAG'] = [40001, 32003, 33004]
    tab['NAME'] = ['CH3CCH', 'CH3OH', 'CH3CH2OH']
    # Q = 10 * T (so the interpolated partition function is easy to predict)
    for ii, tem in enumerate(JPL_TEMPERATURES):
        tab[f'QLOG{ii + 1}'] = np.log10([10 * tem] * 3)
    tab.meta = {'Temperature (K)': JPL_TEMPERATURES}
    return tab


def fake_cdms_species_table():
    tab = Table()
    tab['tag'] = [40001, 66506]
    tab['molecule'] = ['CH3CCH', 'H2C(CN)2']
    for tem in [300., 225., 150., 75., 37.5, 18.75, 9.375]:
        tab[f'lg(Q({tem}))'] = np.log10([10 * tem] * 2)
    tab.meta = {'Temperature (K)': [300., 225., 150., 75., 37.5, 18.75,
                                    9.375]}
    return tab


def fake_catalog_table(tag=40001):
    """A minimal SPCAT-style catalog table like those astroquery returns."""
    tab = Table()
    tab['FREQ'] = np.array([50000., 65000., 105000.]) * u.MHz
    tab['ERR'] = np.array([0.01, 0.01, 0.01]) * u.MHz
    tab['LGINT'] = [-5.0, -4.0, -4.5]
    tab['DR'] = [3, 3, 3]
    tab['ELO'] = np.array([0.0, 5.0, 10.0]) * u.cm**-1
    tab['GUP'] = [3, 5, 7]
    tab['TAG'] = [tag] * 3
    tab['QNFMT'] = [202] * 3
    return tab


@pytest.fixture
def fake_jpl(monkeypatch):
    """Monkeypatch the JPLSpec query tool with fake tables."""
    from astroquery.linelists.jplspec import JPLSpec

    calls = {'get_molecule': [], 'query_lines': []}

    def get_species_table(**kwargs):
        return fake_jpl_species_table()

    def get_molecule(molecule_id, **kwargs):
        calls['get_molecule'].append(molecule_id)
        return fake_catalog_table(tag=molecule_id)

    def query_lines(fmin, fmax, **kwargs):
        calls['query_lines'].append(kwargs)
        return fake_catalog_table()

    monkeypatch.setattr(JPLSpec, 'get_species_table', get_species_table)
    monkeypatch.setattr(JPLSpec, 'get_molecule', get_molecule)
    monkeypatch.setattr(JPLSpec, 'query_lines', query_lines)
    return calls


@pytest.fixture
def fake_cdms(monkeypatch):
    """Monkeypatch the CDMS query tool with fake tables."""
    from astroquery.linelists.cdms import CDMS

    calls = {'get_molecule': [], 'query_lines': []}

    def get_species_table(**kwargs):
        return fake_cdms_species_table()

    def get_molecule(molecule_id, **kwargs):
        calls['get_molecule'].append(molecule_id)
        return fake_catalog_table(tag=molecule_id)

    def query_lines(fmin, fmax, **kwargs):
        calls['query_lines'].append(kwargs)
        return fake_catalog_table()

    monkeypatch.setattr(CDMS, 'get_species_table', get_species_table)
    monkeypatch.setattr(CDMS, 'get_molecule', get_molecule)
    monkeypatch.setattr(CDMS, 'query_lines', query_lines)
    return calls


def test_jpl_name_resolution_and_freq_filter(fake_jpl):
    """
    Regression test for issue #421: the molecule name must be resolved to a
    tag locally and passed to get_molecule, and the returned catalog must be
    restricted to the requested frequency range.
    """
    freqs, aij, deg, EU, partfunc = get_molecular_parameters(
        'CH3CCH', fmin=60 * u.GHz, fmax=110 * u.GHz)

    # name resolved to the correct tag, passed to get_molecule
    assert fake_jpl['get_molecule'] == [40001]
    # JPL default path must not use the (unreliable) query_lines form
    assert fake_jpl['query_lines'] == []

    # the 50 GHz line is out of range; the 65 and 105 GHz lines are kept
    np.testing.assert_allclose(freqs.to(u.MHz).value, [65000., 105000.])
    np.testing.assert_array_equal(deg, [5, 7])
    assert np.all(np.isfinite(aij))
    assert np.all(EU > 0)

    # Q = 10 * T at the tabulated points; 150 K is a node of the interpolant
    assert partfunc(150) == pytest.approx(1500, rel=1e-6)


def test_jpl_tag_selection(fake_jpl):
    freqs, aij, deg, EU, partfunc = get_molecular_parameters(
        None, molecule_tag=32003, fmin=60 * u.GHz, fmax=110 * u.GHz)
    assert fake_jpl['get_molecule'] == [32003]
    assert len(freqs) == 2


def test_jpl_partial_name_match(fake_jpl):
    # 'H3CC' is a substring of only CH3CCH
    freqs, aij, deg, EU, partfunc = get_molecular_parameters(
        'H3CC', fmin=60 * u.GHz, fmax=110 * u.GHz)
    assert fake_jpl['get_molecule'] == [40001]


def test_jpl_regex_name_match(fake_jpl):
    # not a substring of any entry, but a regex matching only CH3CCH
    freqs, aij, deg, EU, partfunc = get_molecular_parameters(
        'CH3..H$', fmin=60 * u.GHz, fmax=110 * u.GHz)
    assert fake_jpl['get_molecule'] == [40001]


def test_jpl_ambiguous_name(fake_jpl):
    # 'CH3' matches all three fake species
    with pytest.raises(ValueError, match='Too many or too few matches'):
        get_molecular_parameters('CH3', fmin=60 * u.GHz, fmax=110 * u.GHz)


def test_jpl_unknown_name(fake_jpl):
    with pytest.raises(ValueError, match='Too many or too few matches'):
        get_molecular_parameters('NOTAMOLECULE',
                                 fmin=60 * u.GHz, fmax=110 * u.GHz)


def test_jpl_no_lines_in_range(fake_jpl):
    with pytest.raises(ValueError, match='No lines found'):
        get_molecular_parameters('CH3CCH', fmin=1 * u.GHz, fmax=2 * u.GHz)


def test_cdms_name_resolution(fake_cdms):
    """
    Regression test for issue #365: names with regex special characters,
    like H2C(CN)2, must resolve; the CDMS default path goes through
    query_lines with a '<zero-padded tag> <name>' search string.
    """
    freqs, aij, deg, EU, partfunc = get_molecular_parameters(
        'H2C(CN)2', catalog='CDMS', fmin=60 * u.GHz, fmax=110 * u.GHz)

    # CDMS default path must use query_lines (astroquery's CDMS
    # catalog-file parser is buggy; see get_molecular_parameters docstring)
    assert fake_cdms['get_molecule'] == []
    assert len(fake_cdms['query_lines']) == 1
    assert fake_cdms['query_lines'][0]['molecule'] == '066506 H2C(CN)2'
    assert fake_cdms['query_lines'][0]['parse_name_locally'] is False

    np.testing.assert_allclose(freqs.to(u.MHz).value, [65000., 105000.])
    assert np.all(np.isfinite(aij))


def test_cdms_freq_filter(fake_cdms):
    freqs, aij, deg, EU, partfunc = get_molecular_parameters(
        'CH3CCH', catalog='CDMS', fmin=60 * u.GHz, fmax=70 * u.GHz)
    assert fake_cdms['query_lines'][0]['molecule'] == '040001 CH3CCH'
    np.testing.assert_allclose(freqs.to(u.MHz).value, [65000.])
    np.testing.assert_array_equal(deg, [5])


def test_cdms_get_molecule_bad_elo(fake_cdms, monkeypatch):
    """
    astroquery's CDMS._parse_cat returns a string-typed ELO column for
    species with letter-coded degeneracies (issue #365); make sure this is
    reported as a clear error instead of a cryptic TypeError.
    """
    from astroquery.linelists.cdms import CDMS

    def get_molecule(molecule_id, **kwargs):
        tab = fake_catalog_table(tag=molecule_id)
        tab.replace_column('ELO', ['0.0', '5.0', '10.0A'])
        return tab

    monkeypatch.setattr(CDMS, 'get_molecule', get_molecule)

    with pytest.raises(ValueError, match='ELO'):
        get_molecular_parameters('H2C(CN)2', catalog='CDMS',
                                 fmin=60 * u.GHz, fmax=110 * u.GHz,
                                 use_get_molecule=True)


@pytest.mark.remote_data
def test_issue421_ch3cch_jpl():
    """https://github.com/pyspeckit/pyspeckit/issues/421"""
    freqs, aij, deg, EU, partfunc = get_molecular_parameters(
        'CH3CCH', fmin=60 * u.GHz, fmax=110 * u.GHz)
    assert len(freqs) > 0
    assert np.all((freqs >= 60 * u.GHz) & (freqs <= 110 * u.GHz))
    # J=4-3, 5-4, 6-5 K-ladders
    assert np.any(np.isclose(freqs.to(u.MHz).value, 102547.9842, atol=1))
    assert np.issubdtype(np.asarray(deg).dtype, np.number)
    assert np.all(np.isfinite(aij))
    assert partfunc(50) > 0


@pytest.mark.remote_data
def test_issue365_h2ccnc2_cdms():
    """https://github.com/pyspeckit/pyspeckit/issues/365"""
    freqs, aij, deg, EU, partfunc = get_molecular_parameters(
        "H2C(CN)2", catalog="CDMS",
        fmin=232567.224454 * u.MHz, fmax=234435.809432 * u.MHz)
    assert len(freqs) > 0
    assert np.issubdtype(np.asarray(deg).dtype, np.number)
    assert np.all(np.isfinite(aij))
    assert np.all(np.isfinite(EU))


@pytest.mark.remote_data
def test_ch3oh_jpl():
    freqs, aij, deg, EU, partfunc = get_molecular_parameters(
        'CH3OH', fmin=60 * u.GHz, fmax=70 * u.GHz)
    assert len(freqs) > 0
    assert np.all((freqs >= 60 * u.GHz) & (freqs <= 70 * u.GHz))
