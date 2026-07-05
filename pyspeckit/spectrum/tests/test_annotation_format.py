"""
Tests for the mathtext annotation formatter (issue #337): plot legends
should show 1.4\\times10^{14} instead of 1.4e+14.
"""
import matplotlib
matplotlib.use('Agg')

import re

import numpy as np
import pytest
from astropy import units as u

from ..annotation_format import format_mathtext_value
from ..classes import Spectrum, units as spunits
from ..models import ammonia, ammonia_constants

# matches computer-style exponents like 'e+14', 'E-23', 'e5'
E_NOTATION = re.compile(r'[eE][-+]?\d')


def assert_no_e_notation(text):
    assert not E_NOTATION.search(text), \
        "found e-notation in {0!r}".format(text)


class TestFormatMathtextValue(object):

    def test_readable_value_stays_fixed_point(self):
        assert format_mathtext_value(3.05) == '3.05'
        assert format_mathtext_value(0.001423) == '0.001423'
        assert format_mathtext_value(99999.0) == '99999'
        assert format_mathtext_value(0) == '0'
        assert format_mathtext_value(-42.5) == '-42.5'

    def test_large_value_uses_times_ten(self):
        out = format_mathtext_value(1.4e14)
        assert out == '1.4\\times10^{14}'
        assert_no_e_notation(out)

    def test_small_value_uses_times_ten(self):
        out = format_mathtext_value(4.89e-23)
        assert out == '4.89\\times10^{-23}'
        assert_no_e_notation(out)

    def test_negative_value_uses_times_ten(self):
        out = format_mathtext_value(-4.95e-19)
        assert out == '-4.95\\times10^{-19}'
        assert_no_e_notation(out)

    def test_shared_exponent_error_form(self):
        out = format_mathtext_value(1.4e14, 5e12)
        assert out == '(1.400\\pm0.050)\\times10^{14}'
        assert_no_e_notation(out)

    def test_shared_exponent_negative_value(self):
        out = format_mathtext_value(-1.4e14, 5e12)
        assert out == '(-1.400\\pm0.050)\\times10^{14}'
        assert_no_e_notation(out)

    def test_readable_value_with_error(self):
        out = format_mathtext_value(4.9, 0.3)
        # two significant figures on the uncertainty
        assert out == '4.90\\pm0.30'
        assert_no_e_notation(out)

    def test_error_none_and_zero(self):
        assert format_mathtext_value(1.4e14, None) == '1.4\\times10^{14}'
        assert format_mathtext_value(1.4e14, 0) == '1.4\\times10^{14}'
        assert format_mathtext_value(5.0, None) == '5'
        assert format_mathtext_value(5.0, 0) == '5'

    def test_value_none(self):
        assert format_mathtext_value(None) == '\\mathrm{None}'

    def test_nan_inf_do_not_crash(self):
        for v in (np.nan, np.inf, -np.inf):
            out = format_mathtext_value(v)
            assert isinstance(out, str)
            out = format_mathtext_value(v, 1.0)
            assert isinstance(out, str)
        out = format_mathtext_value(1.0, np.nan)
        assert isinstance(out, str)
        out = format_mathtext_value(1.0, np.inf)
        assert isinstance(out, str)

    def test_mantissa_rounds_to_ten(self):
        # 9.999e4 -> rounding pushes the mantissa to 10; must renormalize
        out = format_mathtext_value(9.99999e5)
        assert out == '1\\times10^{6}'
        assert_no_e_notation(out)

    @pytest.mark.parametrize('value,error',
                             [(1.4e14, 5e12), (4.89e-23, 1e-25),
                              (-4.95e-19, 3e-21), (3.05e-15, None),
                              (0.0, 0.1), (12.5, 0.02), (np.nan, 1.0),
                              (1e7, 250.), (-99999.0, 1e-7)])
    def test_output_is_valid_mathtext(self, value, error):
        """every output must parse as matplotlib mathtext"""
        from matplotlib.mathtext import MathTextParser
        parser = MathTextParser('agg')
        out = format_mathtext_value(value, error)
        assert_no_e_notation(out)
        parser.parse('$' + out + '$', dpi=72)


def get_legend_texts(sp):
    texts = []
    for legend in sp.plotter.axis.findobj(matplotlib.legend.Legend):
        texts += [t.get_text() for t in legend.get_texts()]
    return texts


def test_gaussian_annotation_renders():
    """
    Fit a Gaussian with a huge amplitude, annotate, and force a draw:
    matplotlib raises on invalid mathtext at draw time.
    """
    import matplotlib.pyplot as plt
    plt.close('all')

    np.random.seed(42)
    x = np.linspace(-50, 50, 201)
    amp = 1.4e14
    data = amp * np.exp(-x**2 / (2 * 5.0**2)) + np.random.randn(201) * amp/100

    sp = Spectrum(xarr=spunits.SpectroscopicAxis(x, unit='km/s'), data=data)
    sp.plotter()
    sp.specfit(fittype='gaussian', guesses=[amp, 0, 5], annotate=True)
    sp.plotter.figure.canvas.draw()

    texts = get_legend_texts(sp)
    assert len(texts) > 0
    joined = "\n".join(texts)
    assert '\\times10^{' in joined
    assert 'e+' not in joined and 'E+' not in joined
    assert_no_e_notation(joined)
    plt.close('all')


def test_ammonia_annotation_renders():
    """
    ammonia_model.annotations is a separate implementation; make sure it
    also produces drawable, e-notation-free labels.
    """
    import matplotlib.pyplot as plt
    plt.close('all')

    velo = u.Quantity(np.linspace(-25, 25, 251), u.km/u.s)
    cfrq = ammonia_constants.freq_dict['oneone'] * u.Hz
    frq = velo.to(u.GHz, u.doppler_radio(cfrq))
    xarr = spunits.SpectroscopicAxis(frq, refX=cfrq)
    data = ammonia.ammonia(xarr, tex=25, trot=25, width=0.5, ntot=13,
                           xoff_v=0.0)
    sp = Spectrum(xarr=xarr, data=data)
    sp.plotter()
    sp.specfit(fittype='ammonia', guesses=[23, 22, 13.1, 1, 0.5, 0],
               fixed=[False, False, False, False, False, True],
               annotate=True)
    sp.plotter.figure.canvas.draw()

    texts = get_legend_texts(sp)
    assert len(texts) > 0
    joined = "\n".join(texts)
    assert 'e+' not in joined and 'E+' not in joined
    assert_no_e_notation(joined)
    plt.close('all')


def test_baseline_annotation_renders():
    """
    The issue #337 screenshot came from the baseline annotation
    (y=+4.89e-23x^2...); make sure it now renders as \\times10^{} form.
    """
    import matplotlib.pyplot as plt
    plt.close('all')

    np.random.seed(0)
    # frequency-like x axis with huge values -> tiny polynomial coefficients
    x = np.linspace(2.15e9, 2.35e9, 200)
    data = 4.89e-23*x**2 - 4.95e-19*x + 3.05e-15 + np.random.randn(200)*1e-6

    sp = Spectrum(xarr=spunits.SpectroscopicAxis(x, unit='Hz'), data=data)
    sp.plotter()
    sp.baseline(order=2, annotate=True, subtract=False)
    sp.plotter.figure.canvas.draw()

    texts = get_legend_texts(sp)
    joined = "\n".join(texts)
    assert '\\times10^{' in joined
    assert 'e+' not in joined and 'E+' not in joined
    assert_no_e_notation(joined)
    plt.close('all')
