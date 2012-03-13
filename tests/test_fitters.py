import pyspeckit
import numpy as np
import pytest

def get_spectrum(offset=25.0, width=2.7, amplitude=3.0, noise=0.0, xlen=101):
    sample_xarr = np.linspace(0,100,xlen)
    sample_data = amplitude * np.exp(-(sample_xarr - offset)**2/(2*width**2)) + np.random.randn(xlen)*noise

    sp = pyspeckit.Spectrum(data=sample_data, xarr=sample_xarr,
            xarr_kwargs={'units':'km/s', 'refX':3e6, 'refX_units':'MHz'})

    return sp

@pytest.mark.parametrize(('use_lmfit','lmfit_type'),
        ((False,'whocares'),
         (True,'leastsq'),
         #(True,'anneal'),
         (True,'lbfgsb')))
def test_fittypes(use_lmfit, lmfit_type):

    sp = get_spectrum()
    if use_lmfit:
        sp.specfit(use_lmfit=use_lmfit, engine=lmfit_type)
    else:
        sp.specfit()

    print sp.specfit.modelpars



