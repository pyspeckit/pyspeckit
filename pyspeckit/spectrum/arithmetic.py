"""
Spectrum Arithmetic

e.g., interpolate one spectrum onto anothers' axes
"""

import classes
import numpy as np

def interp(spec1,spec2):
    """
    Interpolate spec1 onto spec2's axes
    """

    if spec1.xarr.units != spec2.xarr.units:
        spec1.xarr.convert_to_units(spec2.xarr.units)

    newdata = np.interp(spec2.xarr,spec1.xarr,spec1.data)

    if spec1.error is not None:
        newerror = np.interp(spec2.xarr,spec1.xarr,spec1.error) 
    else:
        newerror = None

    newSpec = classes.Spectrum(xarr=spec2.xarr,data=newdata, error=newerror,
            header=spec2.header)

    return newSpec
