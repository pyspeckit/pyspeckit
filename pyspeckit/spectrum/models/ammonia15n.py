"""
========================================
15NH3 inversion transition TROT fitter
========================================

Fitter for the inversion transitions of the 15N-bearing ammonia
isotopologue (15NH3).  This re-uses the NH3 machinery from
`~pyspeckit.spectrum.models.ammonia` with the spectroscopic constants of
15NH3 (`~pyspeckit.spectrum.models.ammonia15n_constants`, based on
calculations by Luca Bizzocchi from a complete reanalysis of the literature
molecular beam data: Kukolich 1967, 1968; Hougen 1972).

.. moduleauthor:: Elena Redaelli <eredaelli@mpe.mpg.de>
"""
from . import ammonia15n_constants
from .ammonia import ammonia, ammonia_model

__all__ = ['ammonia15n', 'ammonia15n_model']


def ammonia15n(xarr, trot=20, tex=None, ntot=12, width=1, xoff_v=0.0,
               fortho=0.0, **kwargs):
    """
    Generate a model 15NH3 spectrum based on input temperatures, column, and
    gaussian parameters.  The returned model will be in Kelvin (brightness
    temperature) units.

    This is `~pyspeckit.spectrum.models.ammonia.ammonia` evaluated with the
    spectroscopic constants of the 15N-bearing isotopologue; see that
    function for the parameter descriptions.  The default column density is
    lower (ntot=12) to reflect the lower abundance of the isotopologue.
    """
    return ammonia(xarr, trot=trot, tex=tex, ntot=ntot, width=width,
                   xoff_v=xoff_v, fortho=fortho,
                   constants=ammonia15n_constants, **kwargs)


class ammonia15n_model(ammonia_model):
    """
    Model and fitter for the 15NH3 inversion transitions, with the same
    parameters as `~pyspeckit.spectrum.models.ammonia.ammonia_model`:
    trot, tex, ntot, width, xoff_v, and fortho.

    Registered as fittype ``ammonia15n``.
    """

    def __init__(self, *args, **kwargs):
        super(ammonia15n_model, self).__init__(*args, **kwargs)
        self.modelfunc = ammonia15n
