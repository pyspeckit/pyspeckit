"""
Regression tests for issue #410: loading an IRAF multispec (echelle)
spectrum, modifying an order, writing it out, and reading it back.

Uses a small synthesized multispec FITS file rather than a real echelle
spectrum so no large data files are needed.
"""
import numpy as np
import pytest
from astropy.io import fits

import pyspeckit
from pyspeckit.wrappers import load_IRAF_multispec

NPIX = 100
# order number, beam, dtype (0=linear), crval, cdelt, npix, z, aplow, aphigh
ORDERS = [(1, 51, 0, 6000.0, 0.5, NPIX, 0.0, 1.0, 2.0),
          (2, 52, 0, 6500.0, 0.4, NPIX, 0.0, 3.0, 4.0)]


def make_multispec_file(filename):
    """Create a minimal IRAF multispec echelle FITS file."""
    rng = np.random.RandomState(42)
    data = rng.rand(len(ORDERS), NPIX).astype('float32')

    hdu = fits.PrimaryHDU(data=data)
    hdr = hdu.header
    hdr['CTYPE1'] = 'MULTISPE'
    hdr['CTYPE2'] = 'MULTISPE'
    hdr['CDELT1'] = 1.0
    hdr['CDELT2'] = 1.0
    hdr['CD1_1'] = 1.0
    hdr['CD2_2'] = 1.0
    hdr['LTM1_1'] = 1.0
    hdr['LTM2_2'] = 1.0
    hdr['WCSDIM'] = 2
    hdr['WAT0_001'] = 'system=multispec'
    hdr['WAT1_001'] = 'wtype=multispec label=Wavelength units=angstroms'
    watspec = ' '.join(['spec{0} = "{0} {1} {2} {3} {4} {5} {6} {7} {8}"'
                        .format(*order) for order in ORDERS])
    watspec = 'wtype=multispec ' + watspec
    # split into 68-character WAT2_xxx cards, as IRAF does
    for ii in range(0, len(watspec), 68):
        hdr['WAT2_%03i' % (ii // 68 + 1)] = watspec[ii:ii + 68]
    hdu.writeto(filename)

    return data


def test_load_write_reload_echelle_order(tmp_path):
    """load_IRAF_multispec -> modify flux -> write -> re-read (issue #410)"""
    msfn = str(tmp_path / 'multispec.fits')
    data = make_multispec_file(msfn)

    orders = load_IRAF_multispec(msfn)
    assert len(orders) == len(ORDERS)

    for sp, order, dd in zip(orders, ORDERS, data):
        crval, cdelt = order[3], order[4]
        np.testing.assert_allclose(np.asarray(sp.xarr.value),
                                   crval + cdelt * np.arange(NPIX))
        np.testing.assert_allclose(np.asarray(sp.data), dd)

    # modify the flux of one order and write it out
    sp = orders[0]
    sp.data = sp.data * 2.0
    outfn = str(tmp_path / 'order0.fits')
    sp.write(outfn, type='fits')

    # writing must not strip the multispec keywords from the in-memory
    # headers (the orders share a single header object)
    assert 'WAT2_001' in sp.header
    assert 'WAT2_001' in orders[1].header

    # the written single-order file must not carry stale multispec keywords
    outhdr = fits.getheader(outfn)
    assert 'WAT0_001' not in outhdr
    assert 'WAT2_001' not in outhdr
    assert 'CD1_1' not in outhdr
    assert outhdr.get('CTYPE2') != 'MULTISPE'

    # ...and must round-trip through the standard 1D reader
    sp2 = pyspeckit.Spectrum(outfn)
    np.testing.assert_allclose(np.asarray(sp2.data), data[0] * 2.0,
                               rtol=1e-6)
    np.testing.assert_allclose(np.asarray(sp2.xarr.value),
                               np.asarray(sp.xarr.value))

    # the error spectrum is stored as the second row
    sp3 = pyspeckit.Spectrum(outfn, errspecnum=1)
    np.testing.assert_allclose(np.asarray(sp3.error),
                               np.asarray(sp.error))


def test_stale_multispec_keywords_error(tmp_path):
    """
    A 2D file with multispec WAT keywords but a non-echelle CTYPE1 (as
    produced by pyspeckit versions that did not strip the stale keywords)
    should raise a comprehensible error, not UnboundLocalError.
    """
    hdu = fits.PrimaryHDU(data=np.zeros((2, NPIX), dtype='float32'))
    hdr = hdu.header
    hdr['CTYPE1'] = 'WAVE'
    hdr['WAT0_001'] = 'system=multispec'
    hdr['WAT1_001'] = 'wtype=multispec label=Wavelength units=angstroms'
    hdr['WAT2_001'] = 'wtype=multispec spec1 = "1 51 0 6000. 0.5 %i 0. 1. 2."' % NPIX
    fn = str(tmp_path / 'stale.fits')
    hdu.writeto(fn)

    with pytest.raises(ValueError, match='MULTISPE'):
        pyspeckit.Spectrum(fn)
