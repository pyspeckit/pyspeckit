"""
Spectrum Arithmetic

e.g., interpolate one spectrum onto anothers' axes

Functions include interp, interp_on_axes, and interpnans.  

 * `interp` will interpolate one spectrum onto another's axes.  
 * `interp_on_axes` will interpolate a spectrum on to a specific X-axis.  
 * `interpnans` will replace NaN values with values interpolated from their
 neighbors.

See the docstrings of these functions for further details.

"""
from __future__ import print_function

import numpy as np

def _interp(x, xp, fp, left=None, right=None):
    """
    Overrides numpy's interp function, which fails to check for
    increasingness....
    """
    if hasattr(x, 'unit'):
        if hasattr(xp, 'unit'):
            x = x.to_value(xp.unit)
            xp = xp.value
        else:
            x = x.value
    if hasattr(xp, 'unit'):
        xp = xp.value
    indices = np.argsort(xp)
    xp = np.array(xp)[indices]
    fp = np.array(fp)[indices]
    return np.interp(x, xp, fp, left, right)

def interp(spec1, spec2, left=0, right=0):
    """
    Interpolate spec1 onto spec2's axes

    Parameters
    ----------
    spec1: pyspeckit.Spectrum
        The Spectrum to interpolate from
    spec2: pyspeckit.Spectrum
        The spectrum to interpolate onto (spec2.xarr will be passed to
        `interp_on_axes`)
    left: float or None
        Passed to np.interp.  Replace values extrapolated on the left
        (negative) side of the output spectrum with this number.
    right: float or None
        Passed to np.interp.  Replace values extrapolated on the right
        (positive) side of the output spectrum with this number.

    Returns
    -------
    spectrum : pyspeckit.Spectrum
        The interpolated spectrum
    """

    return interp_on_axes(spec1, spec2.xarr, left=left, right=right)

def interp_on_axes(spec1, xarr, left=0, right=0):
    """
    Interpolate spec1 onto specified xarr

    Parameters
    ----------
    spec1: pyspeckit.Spectrum
        The Spectrum to interpolate from
    xarr: pyspeckit.xarr
        The X-axis onto which the data will be interpolated
    left: float or None
        Passed to np.interp.  Replace values extrapolated on the left
        (negative) side of the output spectrum with this number.
    right: float or None
        Passed to np.interp.  Replace values extrapolated on the right
        (positive) side of the output spectrum with this number.

    Returns
    -------
    spectrum : pyspeckit.Spectrum
        The interpolated spectrum
    """

    xarr1 = spec1.xarr.as_unit(xarr.unit)

    newdata = _interp(xarr, xarr1, spec1.data, left=left, right=right)

    if spec1.error is not None:
        if xarr1.cdelt() and xarr.cdelt():
            binsizeratio = xarr1.cdelt() / xarr.cdelt()
            # reduce errors by sqrt(binsize) if going from small to large bins
            # else increase errors by similar factor (WARNING!  THEY WILL BE CORRELATED!)
            newerror = _interp(xarr,xarr1,spec1.error) * binsizeratio**0.5
        else:
            newerror = _interp(xarr,xarr1,spec1.error)
    else:
        newerror = None

    newSpec = spec1.copy()
    newSpec.xarr = xarr.copy()
    newSpec.data = newdata
    newSpec.error = newerror

    return newSpec

def interpnans(spec):
    """
    Interpolate over NAN values, replacing them with values interpolated from
    their neighbors using linear interpolation.
    """

    if hasattr(spec.data,'mask'):
        if type(spec.data.mask) is np.ndarray:
            OK = True - spec.data.mask
    if np.any(np.isnan(spec.data) + np.isinf(spec.data)):
        OK = True - (np.isnan(spec.data) + np.isinf(spec.data))

    newdata = _interp(spec.xarr,spec.xarr[OK],spec.data[OK])
    spec.data.mask[:] = False
    spec.data = newdata

    if spec.error is not None:
        newerror = _interp(spec.xarr,spec.xarr[OK],spec.error[OK])
        spec.error = newerror
