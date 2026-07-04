from __future__ import print_function
from .. import astropy_models
import numpy as np
try:
    from astropy.modeling import fitting
    from astropy.modeling import powerlaws
except ImportError:
    pass

# pytest.markif....
def test_powerlaw(scale=5., alpha=2.):
    x = np.linspace(10,100)
    y = scale * (x)**(-alpha)
    plm = powerlaws.PowerLaw1D(amplitude=1, x_0=1, alpha=1)
    # amplitude and x_0 are fully degenerate in PowerLaw1D
    # (amplitude*(x/x_0)**-alpha), and letting both vary can drive the
    # optimizer through x_0 <= 0 where the model is non-finite (which
    # newer astropy fitters reject with NonFiniteValueError)
    plm.x_0.fixed = True
    fitter = fitting.LevMarLSQFitter()
    result = fitter(plm, x, y)
    np.testing.assert_allclose(result.alpha.value, alpha, rtol=1e-6)
    np.testing.assert_allclose(result.amplitude.value, scale, rtol=1e-6)
    return result

if __name__ == "__main__":
    scale,alpha = 5.,2.
    fit_result = test_powerlaw(scale, alpha)
    p = pars = fit_result.parameters

    x = np.linspace(10,100)
    y = scale * (x)**(-alpha)
    yfit = p[0]*(x)**(-p[1])

    import pylab as pl
    pl.plot(x,y)
    pl.plot(x,yfit)
