import pytest
import numpy as np
from astropy import units as u
from ...classes import Spectrum

gg = [5., 2.8, 13.02918494, 0.0855, 14.85, 0.5]


def test_ammonia_parlimits():
    np.random.seed(0)

    sp = Spectrum(xarr=np.linspace(23.68,23.70)*u.GHz,
                  data=np.random.randn(50))

    sp.specfit(fittype='ammonia',
               guesses=[5., 2.8, 13.02918494, 0.0855, 14.85, 0.5],
              )

def test_ammonia_parlimits_fails():
    np.random.seed(0)

    sp = Spectrum(xarr=np.linspace(23.68,23.70)*u.GHz,
                  data=np.random.randn(50))

    with pytest.raises(ValueError) as ex:

        sp.specfit(fittype='ammonia',
                   guesses=[5., 2.8, 13.02918494, 0.0855, 14.85, 0.5],
                   limitedmin=[True, False, False, True, False, True],
                  )

    assert 'no such limit is set' in str(ex.value)
