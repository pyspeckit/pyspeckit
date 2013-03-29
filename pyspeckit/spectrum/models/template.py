from ..interpolation import interp,_interp
import numpy as np

def spectral_template(xarr, template_spectrum, scale, xshift, xshift_units='km/s'):
    """
    Given a template Spectrum (which should be a Spectrum instance),
    scale & shift it
    """
    model = _interp(xarr, template_spectrum.xarr, template_spectrum.data, left=np.nan, right=np.nan)
    pass

