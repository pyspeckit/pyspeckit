import pyspeckit
import numpy as np
import pytest
import timeit
import warnings
warnings.simplefilter("ignore")

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
        sp.specfit(use_lmfit=use_lmfit, engine=lmfit_type, verbose=False)
    else:
        sp.specfit(verbose=False)

    #print sp.specfit.modelpars

sp = get_spectrum()
@pytest.mark.parametrize(('use_lmfit','lmfit_type'),
        ((False,'whocares'),
         (True,'leastsq'),
         #(True,'anneal'),
         (True,'lbfgsb')))
def test_fitspeed(use_lmfit, lmfit_type):
    if use_lmfit:
        sp.specfit(use_lmfit=use_lmfit, engine=lmfit_type, verbose=False)
    else:
        sp.specfit(verbose=False)

    #print sp.specfit.modelpars

lmfit_times = timeit.Timer('test_fitters.test_fittypes(True,"leastsq")',setup='import test_fitters').repeat(5,30)
mpfit_times = timeit.Timer('test_fitters.test_fittypes(False,"leastsq")',setup='import test_fitters').repeat(5,30)

lmfit_nosetup_times = timeit.Timer('test_fitters.test_fitspeed(True,"leastsq")',setup='import test_fitters').repeat(5,30)
mpfit_nosetup_times = timeit.Timer('test_fitters.test_fitspeed(False,"leastsq")',setup='import test_fitters').repeat(5,30)

print "lmfit: %f +/- %f, %f +/- %f" % (np.mean(lmfit_times),np.std(lmfit_times),np.mean(lmfit_nosetup_times),np.std(lmfit_nosetup_times))
print "mpfit: %f +/- %f, %f +/- %f" % (np.mean(mpfit_times),np.std(mpfit_times),np.mean(mpfit_nosetup_times),np.std(mpfit_nosetup_times))
