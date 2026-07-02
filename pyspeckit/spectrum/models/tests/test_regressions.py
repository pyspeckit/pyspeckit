"""
Regression tests for assorted model-layer bugs:

 - hyperfine background fitter parlimited/parlimits misalignment
 - ammonia multi-peak 'tied' strings not being re-indexed
 - _increment_string_number replacing all indices with the first one
 - normalized Gaussian dividing by width**2 instead of width
 - ParinfoList.__setitem__ infinite recursion on existing keys
"""
import numpy as np

from .. import ammonia
from .. import n2hp
from .. import inherited_gaussfitter
from ...parinfo import Parinfo, ParinfoList


def test_hyperfine_background_fitter_parinfo_alignment():
    """
    The background fitters declared 5 parnames but 6-element
    parlimited/parlimits, shifting every limit after 'Tex': tau became
    unbounded, center silently became limited >=0 (forbidding blueshifted
    fits), and width became unbounded.
    """
    for fitter in (n2hp.n2hp_vtau_tbg_fitter,
                   n2hp.n2hp_vtau.background_contsub_fitter):
        parinfo = fitter.make_parinfo()
        limits = {p['parname'].strip('0123456789'): (p['limited'], p['limits'])
                  for p in parinfo}
        assert set(limits.keys()) == {'TBACKGROUND', 'TEX', 'TAU', 'CENTER',
                                      'WIDTH'}
        # Tbackground >= 0
        assert limits['TBACKGROUND'] == ((True, False), (0, 0))
        # Tex >= 1e-5 (Tex=0 is a singularity)
        assert limits['TEX'] == ((True, False), (1e-5, 0))
        # tau >= 0
        assert limits['TAU'] == ((True, False), (0, 0))
        # center must be UNlimited (velocity can be negative)
        assert limits['CENTER'] == ((False, False), (0, 0))
        # width >= 0
        assert limits['WIDTH'] == ((True, False), (0, 0))


def test_increment_string_number():
    """
    _increment_string_number used to replace ALL numbers in the string with
    the (first number + increment): 'p[0]-p[6]' + 7 -> 'p[7]-p[7]'.
    Each index must be incremented independently.
    """
    assert ammonia._increment_string_number('p[0]-p[6]', 7) == 'p[7]-p[13]'
    assert ammonia._increment_string_number('p[6]', 7) == 'p[13]'
    # empty / index-free strings pass through unchanged
    assert ammonia._increment_string_number('', 7) == ''
    assert ammonia._increment_string_number('notied', 7) == 'notied'


def test_ammonia_restricted_tex_multipeak_tied():
    """
    With npeaks=2, the per-peak 'tied' strings used to be replicated
    verbatim, so peak 2's tex was tied to PEAK 1's parameters.
    """
    mod = ammonia.ammonia_model_restricted_tex()
    parinfo = mod.make_parinfo(npeaks=2, params=(20, 20, 0.5, 1.0, 0.0, 0.5, 0)*2)
    assert parinfo.tied == ['', 'p[0]-p[6]', '', '', '', '', '',
                            '', 'p[7]-p[13]', '', '', '', '', '']

    # single-peak behavior must be unchanged
    mod1 = ammonia.ammonia_model_restricted_tex()
    parinfo1 = mod1.make_parinfo()
    assert parinfo1.tied == ['', 'p[0]-p[6]', '', '', '', '', '']


def test_normalized_gaussian_integral():
    """
    The normalized Gaussian used to divide by width**2 instead of width,
    making it wrong by a factor of width.
    """
    # np.trapz was removed in numpy 2.0 in favor of np.trapezoid
    trapezoid = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
    x = np.linspace(-100, 100, 200001)
    for amp, width in [(1.0, 1.0), (3.0, 2.5), (2.0, 0.3)]:
        y = inherited_gaussfitter.gaussian(x, amp, 0.0, width,
                                           normalized=True)
        assert np.isclose(trapezoid(y, x), amp, rtol=1e-6)


def test_parinfolist_setitem():
    """
    ParinfoList.__setitem__ used to recurse infinitely (self[key] = val)
    when replacing an existing item by int or string key.
    """
    pl = ParinfoList([Parinfo({'n': 0, 'value': 1.0, 'parname': 'AMPLITUDE'}),
                      Parinfo({'n': 1, 'value': 2.0, 'parname': 'SHIFT'})])

    # replace by integer index
    pl[0] = Parinfo({'n': 0, 'value': 5.0, 'parname': 'AMPLITUDE0'})
    assert pl.values == [5.0, 2.0]
    assert pl[0]['value'] == 5.0
    assert pl['AMPLITUDE0']['value'] == 5.0

    # replace by string (parameter name) key
    pl['AMPLITUDE0'] = Parinfo({'n': 0, 'value': 7.0,
                                'parname': 'AMPLITUDE0'})
    assert pl.values == [7.0, 2.0]
    assert pl[0]['value'] == 7.0
    assert pl['AMPLITUDE0']['value'] == 7.0

    # out-of-range int assignment still raises
    import pytest
    with pytest.raises(IndexError):
        pl[5] = Parinfo({'n': 5, 'value': 1.0, 'parname': 'NOPE'})
