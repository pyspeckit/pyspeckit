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


def test_cold_ammonia_restricted_tex_multipeak_tied():
    """
    Issue #368: the cold-ammonia restricted-tex tie encodes the nonlinear
    Swift tkin->trot conversion, so its 'tied' string contains numeric
    constants (41.18, 0.6, 15.7).  Multi-peak re-indexing must bump only
    the bracketed parameter indices, not those constants.
    """
    tie0 = ammonia._cold_restricted_tex_tied
    tie1 = ammonia._increment_string_number(tie0, 7)

    mod = ammonia.cold_ammonia_model_restricted_tex()
    parinfo = mod.make_parinfo(npeaks=2,
                               params=(20, 20, 14, 1.0, 0.0, 0.5, 0)*2)
    assert parinfo.tied == ['', tie0, '', '', '', '', '',
                            '', tie1, '', '', '', '', '']

    # the second peak ties tex1 to tkin1 (p[7]) and delta1 (p[13]) ...
    assert 'p[7]' in tie1 and 'p[13]' in tie1
    # ... and the Swift-conversion constants must be untouched
    for constant in ('41.18', '0.6', '15.7'):
        assert constant in tie1


def test_ammonia_restricted_tex_rejects_cold_ammonia_func():
    """
    Issue #368: ammonia_model_restricted_tex(ammonia_func=cold_ammonia)
    was advertised by the docstring but could never work (it crashed with
    TypeError/KeyError deep in the fit); it must now fail up front with a
    pointer to cold_ammonia_model_restricted_tex.
    """
    import pytest
    with pytest.raises(ValueError) as ex:
        ammonia.ammonia_model_restricted_tex(ammonia_func=ammonia.cold_ammonia)
    assert 'cold_ammonia_model_restricted_tex' in str(ex.value)


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


def test_formaldehyde_mm_radex_synthetic_grid():
    """
    formaldehyde_mm_radex indexed its (temperature, column, density) RADEX
    grids with a *list* of slices whose bounds were np.float64; modern numpy
    raises for both (list-of-slices fancy indexing and float slice bounds).
    Feed it a small synthetic grid that is linear in the grid indices, so
    the trilinear interpolation is exact, and compare against the
    directly-evaluated vtau model.
    """
    import astropy.io.fits as pyfits
    from astropy import units as u
    from .. import formaldehyde_mm
    from ...units import SpectroscopicAxis

    ntemp, ncol, ndens = 5, 6, 7
    hdr = pyfits.Header()
    hdr['CRPIX1'], hdr['CRVAL1'], hdr['CDELT1'] = 1, 2.0, 0.5    # log density: 2-5
    hdr['CRPIX2'], hdr['CRVAL2'], hdr['CDELT2'] = 1, 11.0, 0.5   # log column: 11-13.5
    hdr['CRPIX3'], hdr['CRVAL3'], hdr['CDELT3'] = 1, 10.0, 10.0  # temperature: 10-50

    zz, yy, xx = np.indices((ntemp, ncol, ndens), dtype=float)
    taug = 0.01 * (1 + zz + 2 * yy + 3 * xx)
    texg = 5.0 + zz + yy + xx

    xarr = SpectroscopicAxis(np.linspace(218.1, 218.35, 200) * u.GHz)

    # temperature=25 -> gridval3=1.5; column=13.25 -> gridval2=4.5;
    # density=4.25 -> gridval1=4.5 -- all deliberately off-grid so the
    # slice bounds are non-integer floats before conversion.
    spec = formaldehyde_mm.formaldehyde_mm_radex(xarr,
                                                 temperature=25,
                                                 column=13.25,
                                                 density=4.25,
                                                 xoff_v=0.0, width=1.0,
                                                 texgrid=((218., 219., texg),),
                                                 taugrid=((218., 219., taug),),
                                                 hdr=hdr)

    # linear grids => exact trilinear interpolation values
    expected_tau = 0.01 * (1 + 1.5 + 2 * 4.5 + 3 * 4.5)  # 0.25
    expected_tex = 5.0 + 1.5 + 4.5 + 4.5                 # 15.5
    direct = formaldehyde_mm.formaldehyde_mm_vtau(xarr, Tex=expected_tex,
                                                  tau=expected_tau,
                                                  xoff_v=0.0, width=1.0)
    assert np.all(np.isfinite(spec))
    assert spec.max() > 0
    np.testing.assert_allclose(np.asarray(spec), np.asarray(direct),
                               rtol=1e-10, atol=1e-12)


def test_h2co_mm_radex_quantity_comparison():
    """
    h2co_mm_radex compared xarr.as_unit('GHz') (a Quantity) against plain
    floats when masking each line's frequency range, which raises
    UnitConversionError with modern astropy; the comparison must use .value.
    """
    from astropy import units as u
    from .. import h2co_mm
    from ...units import SpectroscopicAxis

    xarr = SpectroscopicAxis(np.linspace(218.1, 218.9, 300) * u.GHz)
    gridbundle = (lambda c, d, t: 12.0,   # Tex303
                  lambda c, d, t: 10.0,   # Tex322
                  lambda c, d, t: 9.0,    # Tex321
                  lambda c, d, t: 0.5,    # tau303
                  lambda c, d, t: 0.2,    # tau322
                  lambda c, d, t: 0.1)    # tau321

    spec = h2co_mm.h2co_mm_radex(xarr, Temperature=25, logColumn=13,
                                 logDensity=4, xoff_v=0.0, width=1.0,
                                 gridbundle=gridbundle)
    assert np.all(np.isfinite(spec))
    assert spec.max() > 0
