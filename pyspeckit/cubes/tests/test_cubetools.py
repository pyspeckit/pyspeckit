from .. import cubes

import numpy as np
from astropy.io import fits
from astropy import wcs

import os

def test_subimage_integ_header():
    # getting a dummy .fits file
    if not os.path.exists('foo.fits'):
        from astropy.utils.data import download_file
        tmp_path = download_file('https://db.tt/oleS9xD6')
        try:
            os.rename(tmp_path, 'foo.fits')
        except OSError:
            # os.rename doesn't like cross-device links
            import shutil
            shutil.move(tmp_path, 'foo.fits')

    cube = fits.getdata('foo.fits')
    header = fits.getheader('foo.fits')

    xcen, ycen = 4.5, 4.5
    xwidth, ywidth = 2.5, 2.5

    # saving results from subimage_integ:
    cutData, cutHead = cubes.subimage_integ(cube, xcen, xwidth, ycen, ywidth,
                                            vrange=(0,9), zunits='pixels',
                                            units='pixels', header=header)

    assert cutHead['CRPIX1'] == 7.0
    assert cutHead['CRPIX2'] == -2.0

    w1 = wcs.WCS(header)
    w2 = wcs.WCS(cutHead)

    # pixel 2,2 in the original image should be pixel 0,0 in the new one
    x1,y1,z1 = w1.wcs_pix2world(2,2,0,0)
    x2,y2 = w2.wcs_pix2world(0,0,0)

    np.testing.assert_almost_equal(x1,x2)
    np.testing.assert_almost_equal(y1,y2)
