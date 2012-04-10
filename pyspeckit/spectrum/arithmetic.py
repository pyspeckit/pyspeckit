"""
Spectrum Arithmetic

e.g., interpolate one spectrum onto anothers' axes
"""

import classes
import numpy as np

def _interp(x, xp, fp, left=None, right=None):
    """
    Overrides numpy's interp function, which fails to check for increasingness....
    """

    if not np.all(np.diff(xp) > 0):
        return np.interp(x, xp[::-1], fp[::-1], left=left, right=right)
    else:
        return np.interp(x, xp, fp, left=left, right=right)

def interp(spec1,spec2):
    """
    Interpolate spec1 onto spec2's axes
    """

    xarr1 = spec1.xarr.as_unit(spec2.xarr.units)

    newdata = _interp(spec2.xarr,xarr1,spec1.data)

    if spec1.error is not None:
        if xarr1.cdelt() and spec2.xarr.cdelt():
            binsizeratio = xarr1.cdelt() / spec2.xarr.cdelt()
            # reduce errors by sqrt(binsize) if going from small to large bins
            # else increase errors by similar factor (WARNING!  THEY WILL BE CORRELATED!)
            newerror = _interp(spec2.xarr,xarr1,spec1.error) * binsizeratio**0.5
        else:
            newerror = _interp(spec2.xarr,xarr1,spec1.error)
    else:
        newerror = None

    newSpec = spec1.copy()
    newSpec.xarr = spec2.xarr.copy()
    newSpec.data = newdata
    newSpec.error = newerror

    return newSpec

def interpnans(spec):
    """
    Interpolate over NAN values, replacing them with their neighbors...
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

