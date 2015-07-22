import warnings
import numpy as np
import pytest

from .. import Spectrum


class TestFitter(object):

    def setup_method(self):

        dx = 0.1
        x = np.arange(-6,6,dx)
        y = 1-np.exp(-x**2 / 2.)
        self.sp = Spectrum(xarr=x, data=y)

    def test_fitter(self):
        self.sp.specfit(fittype='gaussian', guesses=(-1, 0, 0.5))

    def test_invalid_guess(self):
        with pytest.raises(ValueError) as ex:
            self.sp.specfit(fittype='gaussian', guesses=(-1, 0, 0.5),
                       limitedmin=(True,True,True), minvals=(0,0,0))
        assert str(ex.value) == '-1.0 is less than the lower limit 0'

    def test_almost_invalid_guess(self):
        with warnings.catch_warnings(record=True) as w:
        
            self.sp.specfit(fittype='gaussian', guesses=(0-np.spacing(0), 0, 0.5),
                            limitedmin=(True,True,True), minvals=(0,0,0))

        warnings.simplefilter("always")
            
        assert "Guesses have been changed from" in str(w[-1].message)
        assert "is less than the lower limit" in str(w[-2].message)

