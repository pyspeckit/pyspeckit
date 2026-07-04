"""
Regression tests for long-broken legacy I/O, measurement, unit-dict, and
interactive-handler code paths (see the fix-legacy-io-plotting audit).
"""
import warnings

import numpy as np
import pytest

from .. import Spectrum
from ..cosmology import Cosmology
from ..measurements import Measurements
from ..units import SmartCaseNoSpaceDict


def _make_spectrum():
    x = np.linspace(-50, 50, 101)
    y = np.exp(-x**2 / (2 * 5.0**2))
    with warnings.catch_warnings():
        # ignore warning about creating an empty header
        warnings.simplefilter('ignore')
        return Spectrum(xarr=x, data=y, xarrkwargs={'unit': 'km/s'})


def test_cosmology_luminosity_distance():
    """
    scipy.integrate.romberg was removed in scipy 1.15; cosmology.py must use
    quad and return a positive luminosity distance in Mpc (previously it
    crashed on import of romberg, or silently returned 0 without scipy).
    """
    cosm = Cosmology()
    d = cosm.LuminosityDistance(1.0)
    assert np.isfinite(d)
    # a z=1 luminosity distance should be several Gpc for any sane cosmology
    assert 1000 < d < 20000
    # comoving radial distance must be positive too (quad returns (value, err))
    assert cosm.ComovingRadialDistance(0., 1.0) > 0


def test_smart_case_no_space_dict():
    """
    SmartCaseNoSpaceDict.update/fromkeys called __setitem__ with a spurious
    self argument, pop raised an uncaught KeyError instead of forwarding the
    default, setdefault dropped the default, and has_key referenced the
    never-existent dict.has_key__.
    """
    d = SmartCaseNoSpaceDict()

    # update must go through __setitem__ (which also stores the
    # lowercase/no-space alias)
    d.update({'Foo Bar': 1})
    assert d['Foo Bar'] == 1
    assert d['foobar'] == 1
    assert d['FOO BAR'] == 1

    # setdefault: existing key (case-insensitively) returns existing value
    assert d.setdefault('FOO BAR', 99) == 1
    # setdefault: missing key inserts and returns the default
    assert d.setdefault('New Key', 3) == 3
    assert d['newkey'] == 3

    # pop with a default must not raise for missing keys
    assert d.pop('not a key', 'fallback') == 'fallback'
    # pop falls back to the case-normalized key
    assert d.pop('NEWKEY') == 3

    # has_key is case-insensitive and must not crash
    assert d.has_key('foo bar')
    assert not d.has_key('nope')

    # fromkeys builds a new SmartCaseNoSpaceDict with aliases
    d2 = d.fromkeys(['A B'], 5)
    assert isinstance(d2, SmartCaseNoSpaceDict)
    assert d2['A B'] == 5
    assert d2['ab'] == 5


def test_compute_luminosity_single_count():
    """
    compute_luminosity looped over components and added
    compute_flux(pars) * 4 pi d^2 once *per component*, but compute_flux
    already sums over all components, so multi-component lines were
    over-counted by a factor of n_components.
    """
    m = Measurements.__new__(Measurements)
    m.fluxnorm = 1.0
    m.d = 2.0

    # two gaussian components: (amp, center, width) x 2
    pars = [1.0, 5000.0, 2.0, 3.0, 5005.0, 1.0]
    flux = m.compute_flux(pars)
    expected_flux = np.sqrt(2 * np.pi) * (1.0 * 2.0 + 3.0 * 1.0)
    assert np.isclose(flux, expected_flux)

    lum = m.compute_luminosity(pars)
    assert np.isclose(lum, flux * 4. * np.pi * m.d**2)


def test_hdf5_writer_overwrite_false_roundtrip(tmp_path):
    """
    The overwrite=False branch of hdf5_writer never opened the file, so
    writing with overwrite=False crashed with NameError ('f' undefined).
    Also round-trips through the reader, exercising the h5py
    Dataset.value -> [()] fix.
    """
    h5py = pytest.importorskip('h5py')
    from ..readers.hdf5_reader import open_hdf5

    sp = _make_spectrum()
    sp.fileprefix = str(tmp_path / 'spec')

    fn = str(tmp_path / 'spec_out.hdf5')
    sp.write(fn, type='hdf5')
    assert (tmp_path / 'spec_out.hdf5').exists()

    # second write with overwrite=False must pick a new name and actually
    # open the file (previously: NameError)
    sp.write(fn, type='hdf5', overwrite=False)
    newfile = tmp_path / 'spec_fit(1).hdf5'
    assert newfile.exists()

    with h5py.File(str(newfile), 'r') as f:
        np.testing.assert_allclose(f['data'][()], sp.data)
        np.testing.assert_allclose(f['xarr'][()], sp.xarr.value)

    # reader round-trip (no error dataset -> uniform errors)
    data, error, xax, header = open_hdf5(str(newfile))
    np.testing.assert_allclose(data, sp.data)
    np.testing.assert_allclose(error, np.ones_like(sp.data))
    np.testing.assert_allclose(np.asarray(xax), sp.xarr.value)


def _count_event_manager_callbacks(sp):
    n = 0
    cbs = sp.plotter.figure.canvas.callbacks.callbacks
    for eventtype in ('button_press_event', 'key_press_event'):
        for key, val in list(cbs.get(eventtype, {}).items()):
            if hasattr(val, 'func'):
                func = val.func
            elif callable(val):
                func = val()
            else:
                func = None
            if func is not None and 'event_manager' in getattr(func, '__name__', ''):
                n += 1
    return n


def test_interactive_handlers_do_not_accumulate():
    """
    clear_all_connections only recognized callbacks with a .func attribute;
    matplotlib >= 3 stores weakref-like entries, so event_manager callbacks
    were never disconnected and every interactive session stacked another
    set of handlers (interactive selections got out of sync; see issue #346).
    """
    sp = _make_spectrum()
    sp.plotter()

    sp.specfit(interactive=True)
    n_first = _count_event_manager_callbacks(sp)
    assert n_first > 0

    # end the session; all event_manager handlers must be disconnected
    sp.specfit.clear_all_connections()
    assert _count_event_manager_callbacks(sp) == 0

    # a second session must not accumulate handlers
    sp.specfit(interactive=True)
    assert _count_event_manager_callbacks(sp) == n_first
    sp.specfit.clear_all_connections()
    assert _count_event_manager_callbacks(sp) == 0
