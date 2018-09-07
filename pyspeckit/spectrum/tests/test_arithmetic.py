import numpy as np
import warnings

from .. import Spectrum


def test_arithmetic():
    dx = 0.1
    x = np.arange(-6,6,dx)
    y = 1-np.exp(-x**2 / 2.)
    with warnings.catch_warnings():
        # ignore warning about creating an empty header
        warnings.simplefilter('ignore')
        sp = Spectrum(xarr=x, data=y)

    newsp = sp - y

    assert np.all(newsp.data == 0)
