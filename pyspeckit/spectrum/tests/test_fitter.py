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
        self.sp.specfit(fittype='gaussian', guesses=(-1, 0, 0.5),
                        limitedmin=(False,False,True), minpars=(0,0,1e-10))
        assert self.sp.specfit.parinfo.limited == [(False,False),(False,False),(True,False)]
        assert self.sp.specfit.parinfo.limits == [(0,0),(0,0),(1e-10,0)]

    def test_set_pars(self):
        self.sp.specfit(fittype='gaussian', guesses=(-1, 0, 0.5),
                        limitedmin=(True,False,True), minpars=(-5,0,1e-10),
                        fixed=[False,True,False],
                       )
        assert self.sp.specfit.parinfo.limited == [(True,False),(False,False),(True,False)]
        assert self.sp.specfit.parinfo.limits == [(-5,0),(0,0),(1e-10,0)]
        assert self.sp.specfit.parinfo.fixed == [False,True,False]

        self.sp.specfit(fittype='gaussian', guesses=(-1, 0, 0.5),
                        limited=[(True,False),(False,False),(True,False)],
                        limits=[(-5,0),(0,0),(1e-10,0)],
                        fixed=[False,True,False],
                       )
        assert self.sp.specfit.parinfo.limited == [(True,False),(False,False),(True,False)]
        assert self.sp.specfit.parinfo.limits == [(-5,0),(0,0),(1e-10,0)]
        assert self.sp.specfit.parinfo.fixed == [False,True,False]

        self.sp.specfit(fittype='gaussian', guesses=(-1, 0, 0.5),
                        parlimited=[(True,False),(False,False),(True,False)],
                        parlimits=[(-5,0),(0,0),(1e-10,0)],
                        parfixed=[False,True,False],
                       )
        assert self.sp.specfit.parinfo.limited == [(True,False),(False,False),(True,False)]
        assert self.sp.specfit.parinfo.limits == [(-5,0),(0,0),(1e-10,0)]
        assert self.sp.specfit.parinfo.fixed == [False,True,False]

    def test_set_tied(self):
        self.sp.specfit(fittype='gaussian', guesses=(-1, 0, 0.5, -0.1, 0, 0.5),
                        tied=['','','','0.1*p[0]','','']
                       )
        assert self.sp.specfit.parinfo.tied == ['','','','0.1*p[0]','','']

    def test_invalid_guess(self):
        with pytest.raises(ValueError) as ex:
            self.sp.specfit(fittype='gaussian', guesses=(-1, 0, 0.5),
                       limitedmin=(True,True,True), minpars=(0,0,1e-10))
        assert str(ex.value) == '-1.0 is less than the lower limit 0'

    def test_almost_invalid_guess(self):
        with warnings.catch_warnings(record=True) as w:
        
            self.sp.specfit(fittype='gaussian', guesses=(0-np.spacing(0), 0, 0.5),
                            limitedmin=(True,True,True), minpars=(0,0,0))

        warnings.simplefilter("always")
            
        assert "Guesses have been changed from" in str(w[-1].message)
        assert "is less than the lower limit" in str(w[-2].message)

