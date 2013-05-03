from .. import astropy_models
import numpy as np
from astropy.models import fitting

def test_powerlaw(scale=5., alpha=2.):
    x = np.linspace(10,100)
    y = scale * (x)**(-alpha)
    plm = astropy_models.PowerLawModel(1,1,1)
    fitter = fitting.NonLinearLSQFitter(plm)
    result = fitter(x,y)
    print "Result: ",result
    print "plm.params: ",plm
    return plm

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
