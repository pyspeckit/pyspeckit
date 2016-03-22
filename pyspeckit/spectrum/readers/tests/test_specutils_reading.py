try:
    import specutils
    import specutils.io.tests
    SPECUTILS_OK = True
except ImportError:
    SPECUTILS_OK = False

import pytest
import os

from ... import Spectrum

@pytest.mark.skipif(not SPECUTILS_OK)
def test_specutils_aao_reader():

    filename = os.path.join(specutils.io.tests.__path__,
                            'files/AAO.fits')

    sp1d = specutils.io.read_fits.read_fits_spectrum1d(filename)
    sp = Spectrum.from_spectrum1d(sp1d[0])

    assert sp.shape == sp1d[0].flux.shape == (2746,)
    assert sp.xarr.unit.to_string() == 'Angstrom'
    assert all(sp.error==0)
