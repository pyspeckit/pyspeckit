from .. import hyperfine
import numpy as np
from astropy import units as u
from astropy import constants

def test_hyperfine():
    # Regression test for issue 125

    vtau_model = hyperfine.hyperfinemodel(['one'], {'one':0.0},
                                          {'one':100e9},
                                          {'one':1.0},
                                          {'one':1.0})


    xarr = [100]*u.GHz
    Tex = 10.0
    Tbackground = 2.73
    tau = 0.1
    result = vtau_model.hyperfine(xarr, Tex=Tex,
                                  tau=tau, xoff_v=0.0, width=1.0,
                                  Tbackground=Tbackground)

    hoverk = (constants.h.cgs/constants.k_B.cgs).value
    T0 = hoverk * xarr.to(u.Hz).value
    oneminusetotheminustau = (1.0-np.exp(-np.array(tau)))
    term1 = (np.exp(T0/Tex) - 1)**-1
    term2 = (np.exp(T0/Tbackground) - 1)**-1
    correctresult = oneminusetotheminustau*T0*(term1 - term2)

    assert result == correctresult
