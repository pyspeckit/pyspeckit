from .. import cubes, Cube, CubeStack
from ...spectrum.models import n2hp

import numpy as np
from astropy.io import fits
from astropy import wcs

import os

def make_test_cube(shape=(30,9,9), outfile='test.fits', sigma=None, seed=0):
    """
    Generates a simple gaussian cube with noise of
    given shape and writes it as a fits file.
    """
    from astropy.convolution import Gaussian1DKernel, Gaussian2DKernel
    if sigma is None:
        sigma1d, sigma2d = shape[0]/10., np.mean(shape[1:])/5.
    else:
        sigma1d, sigma2d = sigma

    gauss1d = Gaussian1DKernel(stddev=sigma1d, x_size=shape[0])
    gauss2d = Gaussian2DKernel(stddev=sigma2d, x_size=shape[1], 
                                               y_size=shape[2])
    test_cube = gauss1d.array[:,None,None] * gauss2d.array
    test_cube=test_cube/test_cube.max()
    # adding noise:
    np.random.seed(seed)
    noise_cube = (np.random.random(test_cube.shape)-.5)* \
                        np.median(test_cube.std(axis=0))
    test_cube += noise_cube
    true_rms = noise_cube.std()
    
    # making a simple header for the test cube:
    test_hdu = fits.PrimaryHDU(test_cube)
    # the strange cdelt values are a workaround 
    # for what seems to be a bug in wcslib:
    # https://github.com/astropy/astropy/issues/4555
    cdelt1, cdelt2, cdelt3 = -(4e-3+1e-8), 4e-3+1e-8, -0.1
    keylist = {'CTYPE1': 'RA---GLS', 'CTYPE2': 'DEC--GLS', 'CTYPE3': 'VRAD',
               'CDELT1': cdelt1, 'CDELT2': cdelt2, 'CDELT3': cdelt3,
               'CRVAL1': 0, 'CRVAL2': 0, 'CRVAL3': 5, 
               'CRPIX1': 9, 'CRPIX2': 0, 'CRPIX3': 5,
               'CUNIT1': 'deg', 'CUNIT2': 'deg', 'CUNIT3': 'km s-1',
               'BUNIT' : 'K', 'EQUINOX': 2000.0}
    # write out some values used to generate the cube:
    keylist['SIGMA' ] = abs(sigma1d*cdelt3), 'in units of CUNIT3'
    keylist['RMSLVL'] = true_rms
    keylist['SEED'  ] = seed

    test_header = fits.Header()
    test_header.update(keylist)
    test_hdu = fits.PrimaryHDU(data=test_cube, header=test_header)
    test_hdu.writeto(outfile, clobber=True, checksum=True)

def download_test_cube(outfile='test.fits'):
    """
    Downloads a sample fits file from Dropbox (325kB).
    """
    from astropy.utils.data import download_file
    test_cube_url = 'https://db.tt/i0jWA7DU'
    tmp_path = download_file(test_cube_url)
    try:
        os.rename(tmp_path, outfile)
    except OSError:
        # os.rename doesn't like cross-device links
        import shutil
        shutil.move(tmp_path, outfile)

def test_subimage_integ_header(cubefile='test.fits'):
    """
    Checks if the coordinates of the spectral
    cube are drifting away after cropping it.
    """
    # getting a dummy .fits file
    if not os.path.exists(cubefile):
        #download_test_cube(cubefile)
        make_test_cube((100,9,9),cubefile)

    cube = fits.getdata(cubefile)
    header = fits.getheader(cubefile)

    xcen, ycen = 4.5, 4.5
    xwidth, ywidth = 2.5, 2.5

    # saving results from subimage_integ:
    cutData, cutHead = cubes.subimage_integ(cube, xcen, xwidth, ycen, ywidth,
                                            vrange=(0,header['NAXIS3']-1), 
                                            zunits='pixels', units='pixels', 
                                            header=header)

    assert cutHead['CRPIX1'] == 7.0
    assert cutHead['CRPIX2'] == -2.0
    w1 = wcs.WCS(header)
    w2 = wcs.WCS(cutHead)

    # pixel 2,2 in the original image should be pixel 0,0 in the new one
    x1,y1,z1 = w1.wcs_pix2world(2,2,0,0)
    x2,y2 = w2.wcs_pix2world(0,0,0)

    np.testing.assert_almost_equal(x1,x2)
    np.testing.assert_almost_equal(y1,y2)

def do_fiteach(save_cube=None, save_pars=None, show_plot=False):
    """Fits a cube with a gaussian for later use"""
    if save_cube is None:
        save_cube = 'test.fits'
    test_sigma = 10 # in pixel values, each pixel is CDELT3 thick
    make_test_cube((100,10,10), save_cube, 
                      sigma=(test_sigma, 5) )
    # from .. import Cube
    spc = Cube(save_cube)
    guesses = [0.5,0.2,0.8]
    map_rms = np.zeros_like(spc.cube[0])+spc.header['RMSLVL']
    map_snr = spc.cube.std(axis=0)/map_rms
    try:
        import multiprocessing
        ncores = multiprocessing.cpu_count()
    except ImportError:
        ncores = 1
        pass
    spc.fiteach(fittype = 'gaussian',
                guesses = guesses,
                start_from_pixel = (5,5),
                multicore = ncores,
                blank_value = np.nan,
                verbose_level = 3,
                errmap = map_rms,
                signal_cut = 5)
    if show_plot:
        spc.mapplot()
    if save_pars:
        spc.write_fit(save_pars, clobber=True)
    return spc

def test_fiteach(save_cube=None, save_pars=None, show_plot=False):
    """
    A simple test on Cube.fiteach() checking
    that for a noise with set seed the fraction
    of line width values within errorbars is
    remaning constant.
    """
    spc = do_fiteach(save_cube, save_pars, show_plot)

    # checking the fit
    map_seed = spc.header['SEED']
    map_sigma_post = spc.parcube[2]
    map_sigma_true = np.zeros_like(map_sigma_post) + spc.header['SIGMA']
    map_in_bounds = np.abs(map_sigma_true-map_sigma_post) < spc.errcube[2]
    err_frac = map_in_bounds[~map_in_bounds].size / float(map_sigma_post.size)

    assert map_seed == 0
    assert err_frac == 0.39

def test_get_modelcube(cubefile=None, parfile=None, sigma_threshold=5):
    """
    Tests get_modelcube() method for Cube and CubeStack classes.
    If either cubefile or parfile isn't set, fill generate and
    fit a sample cube through do_fiteach().

    Computes the residual cube and collapses it into standard
    deviation of the residual map. Checks that the number of the
    residual pixels three sigma doesn't change for a fixed noise.
    """
    if cubefile is None or parfile is None:
        cubefile = 'test.fits'
        parfile  = 'test_pars.fits'
        sp_cube = do_fiteach(save_cube=cubefile, save_pars=parfile)
    else:
        sp_cube  = Cube(cubefile)
    map_seed = sp_cube.header['SEED']
    map_rms = sp_cube.header['RMSLVL']
    sp_cube.xarr.velocity_convention = 'radio'
    sp_stack = CubeStack([sp_cube])
    sp_stack._modelcube = None
    # assuming one gaussian component
    for spc in [sp_cube, sp_stack]:
        spc.load_model_fit(parfile, npars=3)
        spc.get_modelcube()
        resid_cube = spc.cube - spc._modelcube
    above3sig = (resid_cube.std(axis=0) > map_rms*3).flatten()

    assert map_seed == 0
    assert above3sig[above3sig].size == 77

def test_get_modelcube_badpar(cubefile=None, parfile=None, sigma_threshold=5):
    """
    Test loading a model cube that has at least one invalid parameter.
    Regression test for #163

    This is essentially only testing that get_modelcube works in the presence
    of invalid fit parameters
    """

    if cubefile is None or parfile is None:
        cubefile = 'test.fits'
        fh = fits.open('test_pars.fits')
        fh[0].data[1,0,0] *= -1 # set the width to be negative
        fh.writeto('test_pars_bad.fits', clobber=True)
        parfile  = 'test_pars_bad.fits'
        sp_cube = do_fiteach(save_cube=cubefile, save_pars=parfile)
    else:
        sp_cube  = Cube(cubefile)
    map_seed = sp_cube.header['SEED']
    map_rms = sp_cube.header['RMSLVL']
    sp_cube.xarr.velocity_convention = 'radio'
    sp_stack = CubeStack([sp_cube])
    sp_stack._modelcube = None
    # assuming one gaussian component
    for spc in [sp_cube, sp_stack]:
        spc.load_model_fit(parfile, npars=3, _temp_fit_loc=(0,0))
        spc.get_modelcube()
        resid_cube = spc.cube - spc._modelcube

def test_registry_inheritance(cubefile='test.fits'):
    """
    Regression test for #166
    """
    # getting a dummy .fits file
    if not os.path.exists(cubefile):
        #download_test_cube(cubefile)
        make_test_cube((100,9,9),cubefile)

    spc = Cube(cubefile)

    spc.xarr.velocity_convention = 'radio'
    # spc.Registry.add_fitter('n2hp_vtau', n2hp.n2hp_vtau_fitter, 4)

    sp = spc.get_spectrum(3,3)
    sp.Registry.add_fitter('n2hp_vtau', n2hp.n2hp_vtau_fitter, 4)
    assert 'n2hp_vtau' in sp.Registry.multifitters
    assert 'n2hp_vtau' in sp.Registry.npars
    sp.specfit(fittype='n2hp_vtau', guesses=[1,2,3,4])
