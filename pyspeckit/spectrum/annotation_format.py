"""
Helpers for formatting numeric values (and their uncertainties) as
matplotlib mathtext for use in plot annotations (fit-parameter legends,
baseline labels, etc.).

The goal is to avoid computer-style exponential notation like ``1.4e+14``
or ``4.89e-23`` in plot legends and instead render such values as
:math:`1.4\\times10^{14}`, which is both prettier and less confusing
(see issue #337).

All returned strings are valid matplotlib mathtext *without* the
enclosing dollar signs; callers are expected to wrap them in ``$...$``
themselves.
"""
from __future__ import print_function

import math

import numpy as np

__all__ = ['format_mathtext_value']

# Values whose absolute value falls in [READABLE_MIN, READABLE_MAX) are
# shown in plain fixed-point notation; everything else uses
# mantissa \times 10^{exponent} notation.
READABLE_MIN = 1e-3
READABLE_MAX = 1e5

# never show more than this many digits after the decimal point
MAX_DECIMALS = 8


def _is_readable(value):
    """Can ``value`` be displayed in fixed-point notation legibly?"""
    return value == 0 or (READABLE_MIN <= abs(value) < READABLE_MAX)


def _exponent(value):
    """Base-10 exponent of ``value`` (value must be finite & nonzero)."""
    return int(math.floor(math.log10(abs(value))))


def _decimals_for(reference):
    """
    Number of digits to display after the decimal point such that
    ``reference`` (usually an uncertainty) is shown with two significant
    figures.  Mirrors the historical ``Decimal(...).quantize('%0.2g')``
    behavior.
    """
    if reference == 0 or not np.isfinite(reference):
        return 2
    return min(MAX_DECIMALS, max(0, 1 - _exponent(reference)))


def _format_sci(value, sigfig=4):
    """Format a finite, nonzero value as ``m\\times10^{e}`` mathtext."""
    exponent = _exponent(value)
    mantissa = value / 10.**exponent
    mstr = '{0:.{1}g}'.format(mantissa, sigfig)
    # rounding can push the mantissa to +/-10; renormalize if so
    if abs(float(mstr)) >= 10:
        exponent += 1
        mantissa = value / 10.**exponent
        mstr = '{0:.{1}g}'.format(mantissa, sigfig)
    return '{0}\\times10^{{{1}}}'.format(mstr, exponent)


def _format_single(value):
    """Format a single value (no uncertainty) as mathtext."""
    if value is None:
        return '\\mathrm{None}'
    value = float(value)
    if np.isnan(value):
        return '\\mathrm{NaN}'
    if np.isinf(value):
        return '\\infty' if value > 0 else '-\\infty'
    if _is_readable(value):
        return '{0:.6g}'.format(value)
    return _format_sci(value)


def format_mathtext_value(value, error=None):
    r"""
    Format a value and optional symmetric uncertainty as matplotlib
    mathtext (without enclosing ``$``).

    Values in the "readable" range (``1e-3 <= |v| < 1e5``, or exactly 0)
    are rendered in fixed-point notation; other values are rendered as
    ``mantissa\times10^{exponent}``.  When both a value and an error are
    given and the value requires exponential notation, the shared-exponent
    form ``(m\pm e)\times10^{exp}`` is used, with the error scaled to the
    value's exponent.

    Parameters
    ----------
    value : float or None
        The value to format.
    error : float or None
        The symmetric 1-sigma uncertainty.  ``None`` or ``0`` means "do
        not display an uncertainty".

    Returns
    -------
    str
        A valid matplotlib mathtext fragment (no ``$`` delimiters).

    Examples
    --------
    >>> print(format_mathtext_value(1.4e14, 5e12))
    (1.400\pm0.050)\times10^{14}
    >>> print(format_mathtext_value(0.001423))
    0.001423
    >>> print(format_mathtext_value(1.4e-3, 5e-5))
    0.001400\pm0.000050
    """
    if value is None or (np.isscalar(value) and not np.isfinite(value)):
        # NaN / inf / None values: show them literally; append the error
        # if one is meaningful
        valstr = _format_single(value)
        if error in (None, 0):
            return valstr
        return '{0}\\pm{1}'.format(valstr, _format_single(abs(error)))

    value = float(value)

    if error in (None, 0):
        return _format_single(value)

    error = abs(float(error))

    if not np.isfinite(error):
        return '{0}\\pm{1}'.format(_format_single(value),
                                   _format_single(error))

    if _is_readable(value):
        # the error may be displayed in fixed-point notation even if it is
        # below READABLE_MIN, as long as it doesn't require more than
        # MAX_DECIMALS digits after the decimal point
        error_fixed_ok = (10.**(1 - MAX_DECIMALS) <= error < READABLE_MAX)
        if error_fixed_ok:
            # fixed-point for both; the number of displayed digits is set
            # by the smaller of |value|, error (historical behavior)
            ref = min(abs(value), error) or max(abs(value), error)
            dv = _decimals_for(ref)
            de = _decimals_for(error)
            return '{0:.{1}f}\\pm{2:.{3}f}'.format(value, dv, error, de)
        else:
            # error requires exponential notation but value does not
            return '{0}\\pm{1}'.format(_format_single(value),
                                       _format_single(error))

    # value requires exponential notation: use the shared-exponent form
    # (m \pm e) x 10^{exp} with the error scaled to the value's exponent
    exponent = _exponent(value)
    scaled_value = value / 10.**exponent
    scaled_error = error / 10.**exponent
    decimals = _decimals_for(scaled_error)
    # rounding may push the mantissa to +/-10; renormalize if so
    if abs(round(scaled_value, decimals)) >= 10:
        exponent += 1
        scaled_value = value / 10.**exponent
        scaled_error = error / 10.**exponent
        decimals = _decimals_for(scaled_error)
    return '({0:.{2}f}\\pm{1:.{2}f})\\times10^{{{3}}}'.format(
        scaled_value, scaled_error, decimals, exponent)
