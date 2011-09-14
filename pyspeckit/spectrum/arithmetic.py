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

    if spec1.xarr.units != spec2.xarr.units:
        spec1.xarr.convert_to_unit(spec2.xarr.units)

    newdata = _interp(spec2.xarr,spec1.xarr,spec1.data)

    if spec1.error is not None:
        newerror = _interp(spec2.xarr,spec1.xarr,spec1.error) 
    else:
        newerror = None

    newSpec = classes.Spectrum(xarr=spec2.xarr,data=newdata, error=newerror,
            header=spec1.header) # inherit the old spectrum's header info
                                 # but nearly all other properties are changed

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

