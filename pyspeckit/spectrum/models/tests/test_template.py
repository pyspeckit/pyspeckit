"""
Tests for template fitter
"""

from pyspeckit import spectrum, models
from pyspeckit.spectrum.classes import Spectrum
import numpy as np

def test_template():
    xarr = spectrum.units.SpectroscopicAxis(np.linspace(-5,5,100),
                                                      unit='km/s')
    gauss = np.exp(-xarr.value**2/(2.*1.**2))
    gauss2 = np.exp(-(xarr.value-1)**2/(2.*1.**2))
    np.random.seed(0)
    noise = np.random.randn(xarr.size) / 100.

    sp = Spectrum(xarr=xarr, data=gauss2+noise, header={})
    template = Spectrum(xarr=xarr, data=gauss, header={})

    template_fitter = models.template_fitter(template,
                                                       xshift_units='km/s')
    sp.Registry.add_fitter('template', template_fitter, 2)
    sp.specfit(fittype='template', guesses=[1,0])

    np.testing.assert_almost_equal(sp.specfit.parinfo.SCALE0.value, 1.00, 2)
    np.testing.assert_almost_equal(sp.specfit.parinfo.SHIFT0.value, 1.00, 2)

    return sp

def test_template_withcont():
    xarr = spectrum.units.SpectroscopicAxis(np.linspace(-5,5,100),
                                                      unit='km/s')

    scale = 0.5
    shift = 1.1

    gauss = np.exp(-xarr.value**2/(2.*1.**2))
    cont_ = np.linspace(0.5,1,100)
    cont = np.linspace(-5,5,100)*0.05+0.75
    np.testing.assert_array_almost_equal(cont_, cont)

    gauss2 = np.exp(-(xarr.value-shift)**2/(2.*1.**2))
    # 0.5 = m*-5 + b
    # 1.0 = m*5 + b
    # 0.5 = m*-5 + 1 - m*5
    # -0.5 = -10*m
    # m = 0.05
    # b = 0.5 + 5*0.05 = 0.75
    contshift = np.linspace(-5-shift,5-shift,100)*0.05+0.75
    np.random.seed(0)
    noise = np.random.randn(xarr.size) / 100.

    sp = Spectrum(xarr=xarr, data=(gauss2+contshift+noise)*scale, header={})
    template = Spectrum(xarr=xarr, data=gauss+cont, header={})

    template_fitter = models.template_fitter(template,
                                                       xshift_units='km/s')
    sp.Registry.add_fitter('template', template_fitter, 2)
    sp.specfit(fittype='template', guesses=[scale,shift])

    np.testing.assert_almost_equal(sp.specfit.parinfo.SCALE0.value, scale, 2)
    np.testing.assert_almost_equal(sp.specfit.parinfo.SHIFT0.value, shift, 2)

    return sp,template_fitter


if __name__ == "__main__":

    sp = test_template()
    sp.plotter()
    sp.specfit.plot_fit()


    sp2,template2 = test_template_withcont()
    sp2.plotter()
    sp2.specfit.plot_fit()
    #mod = template2.modelfunc(sp2.xarr, 0.5, 1.1)
    #templ = template2.modelfunc(sp2.xarr, 1, 0)
    #import pylab as pl
    #sp2.plotter.axis.plot(sp2.xarr, mod, 'b')
    #sp2.plotter.axis.plot(sp2.xarr, 0.5*templ, 'g')

    #cont = np.linspace(-5,5,100)*0.05+0.75
    #shift=1.1
    #contshift = np.linspace(-5+shift,5+shift,100)*0.05+0.75
    #sp2.plotter.axis.plot(sp2.xarr, 0.5*contshift, 'b')
    #sp2.plotter.axis.plot(sp2.xarr, 0.5*cont, 'g')

    #pl.draw()
    #pl.show()

