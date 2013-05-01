from .. import astropy_models
import numpy as np
from astropy.models import fitting

def test_powerlaw(scale=1., alpha=1., characteristic=10.):
    x = np.linspace(10,100)
    y = scale * (x/characteristic)**(-alpha)
    plm = astropy_models.PowerLawModel(0,0,0)
    fitter = fitting.NonLinearLSQFitter(plm)
    result = fitter(x,y)
    print "Result: ",result
    print "plm.params: ",plm

if __name__ == "__main__":
    test_powerlaw()
