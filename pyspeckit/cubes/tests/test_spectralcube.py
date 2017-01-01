import os
import numpy as np

from .. import Cube
from .test_cubetools import make_test_cube


def test_get_spectrum(cubefile='smalltest.fits'):
    # getting a dummy .fits file
    if not os.path.exists(cubefile):
        #download_test_cube(cubefile)
        make_test_cube((10,3,3),cubefile)

    pcube = Cube(cubefile)

    sp = pcube.get_spectrum(1,1)

    assert len(sp) == 10
