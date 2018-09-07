"""
This test is skipped until specutils is restored to working status
"""


# try:
#     import specutils
#     import specutils.io.tests
#     SPECUTILS_OK = True
# except ImportError:
#     SPECUTILS_OK = False
# 
# import pytest
# import os
# 
# from ... import Spectrum, Spectra
# 
# @pytest.mark.skipif("not SPECUTILS_OK")
# def test_specutils_aao_reader_single():
# 
#     filename = os.path.join(specutils.io.tests.__path__[0],
#                             'files/AAO.fits')
# 
#     sp1d = specutils.io.read_fits.read_fits_spectrum1d(filename)
#     sp = Spectrum.from_spectrum1d(sp1d[0])
# 
#     assert sp.shape == sp1d[0].flux.shape == (2746,)
#     assert sp.xarr.unit.to_string() == 'Angstrom'
#     assert all(sp.error==0)
# 
# @pytest.mark.skipif("not SPECUTILS_OK")
# def test_specutils_aao_reader_multiple():
# 
#     filename = os.path.join(specutils.io.tests.__path__[0],
#                             'files/AAO.fits')
# 
#     sp1d = specutils.io.read_fits.read_fits_spectrum1d(filename)
#     sp = Spectra.from_spectrum1d_list(sp1d)
# 
#     assert sp.shape[0] == sum([x.flux.shape[0] for x in sp1d]) == 140046
#     assert sp.xarr.unit.to_string() == 'Angstrom'
#     assert all(sp.error==0)
# 
# @pytest.mark.skipif("not SPECUTILS_OK")
# def test_specutils_aao_reader_dontallowmismatchdiffs():
# 
#     filename = os.path.join(specutils.io.tests.__path__[0],
#                             'files/AAO.fits')
# 
#     sp1d = specutils.io.read_fits.read_fits_spectrum1d(filename)
#     sp = Spectra.from_spectrum1d_list(sp1d)
# 
#     # TODO: make this an independent test
#     # this is testing that you can't do arithmetic on axes that are different
#     with pytest.raises(ValueError) as ex:
#         sp[0]-sp[1]
#     assert str(ex.value) == "X-axes do not match."
