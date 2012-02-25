"""
A simple power-law extinction model

Returns multiplicative factors...
"""
from pyspeckit.spectrum import units

def extinction_function(wavelength, tau0=0, alpha=0, lambda0=2.5, lambda_units='microns'):
    """
    Returns a model with modelfunc replaced...
    """
    if hasattr(wavelength,'as_unit'):
        wavelength = wavelength.as_unit(lambda_units)

    tau_of_lambda = tau0 * ( (wavelength / lambda0)**alpha - 1.0)

    return np.exp(-tau_of_lambda)
