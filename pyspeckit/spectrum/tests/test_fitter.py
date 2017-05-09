import warnings
import numpy as np
import pytest

from .. import Spectrum


class TestFitter(object):

    def setup_method(self):

        dx = 0.1
        x = np.arange(-6,6,dx)
        y = -np.exp(-x**2 / 2.)
        with warnings.catch_warnings():
            # ignore warning about creating an empty header
            warnings.simplefilter('ignore')
            self.sp = Spectrum(xarr=x, data=y)

        y_multi = np.exp(-(x-3)**2/2.) - np.exp(-(x+3)**2/2.)
        with warnings.catch_warnings():
            # ignore warning about creating an empty header
            warnings.simplefilter('ignore')
            self.sp_multi = Spectrum(xarr=x, data=y_multi)

        y_multi_tiny = 1e-12 * y_multi
        with warnings.catch_warnings():
            # ignore warning about creating an empty header
            warnings.simplefilter('ignore')
            self.sp_multi_tiny = Spectrum(xarr=x, data=y_multi_tiny)


    def test_init(self):
        assert not all(self.sp.specfit.errspec == 0)

    def test_copy(self):
        """
        Regression test for #166: make sure that a Spectrum's Registry is the
        same as its fitter's registry
        """
        spcopy = self.sp.copy()
        assert self.sp.Registry is self.sp.specfit.Registry
        assert spcopy.Registry is spcopy.specfit.Registry

    def test_fitter(self, use_lmfit=False):
        self.sp.specfit(fittype='gaussian', guesses=(-1, 0, 0.5),
                        limitedmin=(False,False,True), minpars=(0,0,1e-10),
                        use_lmfit=use_lmfit)
        assert self.sp.specfit.parinfo.limited == [(False,False),(False,False),(True,False)]
        assert self.sp.specfit.parinfo.limits == [(0,0),(0,0),(1e-10,0)]
        np.testing.assert_almost_equal(self.sp.specfit.parinfo[0].value, -1)
        np.testing.assert_almost_equal(self.sp.specfit.parinfo[1].value,  0)
        np.testing.assert_almost_equal(self.sp.specfit.parinfo[2].value,  1)

        # Do it again to make sure there are no behavioral changes the second time around
        self.sp.specfit(fittype='gaussian', guesses=(-1, 0, 0.5),
                        limitedmin=(False,False,True), minpars=(0,0,1e-10),
                        use_lmfit=use_lmfit)
        assert self.sp.specfit.parinfo.limited == [(False,False),(False,False),(True,False)]
        assert self.sp.specfit.parinfo.limits == [(0,0),(0,0),(1e-10,0)]
        np.testing.assert_almost_equal(self.sp.specfit.parinfo[0].value, -1)
        np.testing.assert_almost_equal(self.sp.specfit.parinfo[1].value,  0)
        np.testing.assert_almost_equal(self.sp.specfit.parinfo[2].value,  1)

    def test_lmfitter(self):
        """
        Regression test for #223: make sure lmfitter returns proper pars
        """
        self.test_fitter(use_lmfit=True)

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
            warnings.simplefilter('default')

            self.sp.specfit(fittype='gaussian', guesses=(0-np.spacing(0), 0, 0.5),
                            limitedmin=(True,True,True), minpars=(0,0,0))

        assert "Guesses have been changed from" in str(w[-1].message)
        assert "is less than the lower limit" in str(w[-2].message)

    def test_multipeak(self):
        self.sp_multi.specfit(fittype='gaussian', guesses=(0.9, 3, 0.5, -0.9, -3,
                                                           0.5),
                              limitedmin=(False,False,True)*2,
                              minpars=(0,0,1e-10)*2)
        assert self.sp_multi.specfit.parinfo.limited == [(False,False),(False,False),(True,False)]*2
        assert self.sp_multi.specfit.parinfo.limits == [(0,0),(0,0),(1e-10,0)]*2

        np.testing.assert_almost_equal(self.sp_multi.specfit.parinfo[0].value,  1)
        np.testing.assert_almost_equal(self.sp_multi.specfit.parinfo[1].value,  3)
        np.testing.assert_almost_equal(self.sp_multi.specfit.parinfo[2].value,  1)

        np.testing.assert_almost_equal(self.sp_multi.specfit.parinfo[3].value, -1)
        np.testing.assert_almost_equal(self.sp_multi.specfit.parinfo[4].value, -3)
        np.testing.assert_almost_equal(self.sp_multi.specfit.parinfo[5].value,  1)

        # Do it again to make sure there are no behavioral changes the second time around
        self.sp_multi.specfit(fittype='gaussian', guesses=(0.9, 3, 0.5, -0.9, -3,
                                                           0.5),
                              limitedmin=(False,False,True)*2,
                              minpars=(0,0,1e-10)*2)
        assert self.sp_multi.specfit.parinfo.limited == [(False,False),(False,False),(True,False)]*2
        assert self.sp_multi.specfit.parinfo.limits == [(0,0),(0,0),(1e-10,0)]*2

        np.testing.assert_almost_equal(self.sp_multi.specfit.parinfo[0].value,  1)
        np.testing.assert_almost_equal(self.sp_multi.specfit.parinfo[1].value,  3)
        np.testing.assert_almost_equal(self.sp_multi.specfit.parinfo[2].value,  1)

        np.testing.assert_almost_equal(self.sp_multi.specfit.parinfo[3].value, -1)
        np.testing.assert_almost_equal(self.sp_multi.specfit.parinfo[4].value, -3)
        np.testing.assert_almost_equal(self.sp_multi.specfit.parinfo[5].value,  1)

    def test_multipeak_tiny(self):
        self.sp_multi_tiny.specfit(fittype='gaussian', guesses=(0.9e-12, 3, 0.5, -0.9e-12, -3,
                                                           0.5),
                                   limitedmin=(False,False,True)*2,
                                   minpars=(0,0,1e-10)*2)
        assert self.sp_multi_tiny.specfit.parinfo.limited == [(False,False),(False,False),(True,False)]*2
        assert self.sp_multi_tiny.specfit.parinfo.limits == [(0,0),(0,0),(1e-10,0)]*2

        np.testing.assert_almost_equal(self.sp_multi_tiny.specfit.parinfo[0].value,  1e-12)
        np.testing.assert_almost_equal(self.sp_multi_tiny.specfit.parinfo[1].value,  3)
        np.testing.assert_almost_equal(self.sp_multi_tiny.specfit.parinfo[2].value,  1)

        np.testing.assert_almost_equal(self.sp_multi_tiny.specfit.parinfo[3].value, -1e-12)
        np.testing.assert_almost_equal(self.sp_multi_tiny.specfit.parinfo[4].value, -3)
        np.testing.assert_almost_equal(self.sp_multi_tiny.specfit.parinfo[5].value,  1)

        # Do it again to make sure there are no behavioral changes the second time around
        self.sp_multi_tiny.specfit(fittype='gaussian', guesses=(0.9e-12, 3, 0.5, -0.9e-12, -3,
                                                           0.5),
                                   limitedmin=(False,False,True)*2,
                                   minpars=(0,0,1e-10)*2)
        assert self.sp_multi_tiny.specfit.parinfo.limited == [(False,False),(False,False),(True,False)]*2
        assert self.sp_multi_tiny.specfit.parinfo.limits == [(0,0),(0,0),(1e-10,0)]*2

        np.testing.assert_almost_equal(self.sp_multi_tiny.specfit.parinfo[0].value,  1e-12)
        np.testing.assert_almost_equal(self.sp_multi_tiny.specfit.parinfo[1].value,  3)
        np.testing.assert_almost_equal(self.sp_multi_tiny.specfit.parinfo[2].value,  1)

        np.testing.assert_almost_equal(self.sp_multi_tiny.specfit.parinfo[3].value, -1e-12)
        np.testing.assert_almost_equal(self.sp_multi_tiny.specfit.parinfo[4].value, -3)
        np.testing.assert_almost_equal(self.sp_multi_tiny.specfit.parinfo[5].value,  1)
