from ..interpolation import interp,_interp
import numpy as np

def spectral_template(xarr, template_spectrum, scale, xshift, xshift_units='km/s'):
    """
    Given a template Spectrum (which should be a Spectrum instance),
    scale & shift it
    """
    shift = xarr.x_to_units(xshift, xshift_units)
    model = _interp(xarr, template_spectrum.xarr.as_unit(xarr.units)+shift,
            template_spectrum.data, left=np.nan, right=np.nan) * scale
    
    return model

