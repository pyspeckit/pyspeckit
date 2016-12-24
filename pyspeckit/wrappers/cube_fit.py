"""
Cube Fitting
============

Complicated code for fitting of a whole data cube, pixel-by-pixel
"""
from __future__ import print_function
import pyspeckit
from astropy.io import fits
import astropy
import numpy as np
import os

def cube_fit(cubefilename, outfilename, errfilename=None, scale_keyword=None,
        vheight=False, verbose=False, signal_cut=3, verbose_level=2,
        overwrite=True, **kwargs):
    """
    Light-weight wrapper for cube fitting

    Takes a cube and error map (error will be computed naively if not given)
    and computes moments then fits for each spectrum in the cube.  It then
    saves the fitted parameters to a reasonably descriptive output file whose
    header will look like ::

        PLANE1  = 'amplitude'
        PLANE2  = 'velocity'
        PLANE3  = 'sigma'
        PLANE4  = 'err_amplitude'
        PLANE5  = 'err_velocity'
        PLANE6  = 'err_sigma'
        PLANE7  = 'integral'
        PLANE8  = 'integral_error'
        CDELT3  = 1
        CTYPE3  = 'FITPAR'
        CRVAL3  = 0
        CRPIX3  = 1

    Parameters
    ----------
    errfilename: [ None | string name of .fits file ]
        A two-dimensional error map to use for computing signal-to-noise cuts
    scale_keyword: [ None | Char ]
        Keyword to pass to the data cube loader - multiplies cube by the number
        indexed by this header kwarg if it exists.  e.g., if your cube is in
        T_A units and you want T_A*
    vheight: [ bool ]
        Is there a background to be fit?  Used in moment computation
    verbose: [ bool ]
    verbose_level: [ int ]
        How loud will the fitting procedure be?  Passed to momenteach and fiteach
    signal_cut: [ float ]
        Signal-to-Noise ratio minimum.  Spectra with a peak below this S/N ratio
        will not be fit and will be left blank in the output fit parameter cube
    overwrite: [ bool ]
        Overwrite parameter .fits cube if it exists?

    `kwargs` are passed to :class:`pyspeckit.Spectrum.specfit`
    """

    # Load the spectrum
    sp = pyspeckit.Cube(cubefilename,scale_keyword=scale_keyword)
    if os.path.exists(errfilename):
        errmap = fits.getdata(errfilename)
    else:
        # very simple error calculation... biased and bad, needs improvement
        # try using something from the cubes package
        errmap = sp.cube.std(axis=0)

    # Run the fitter
    sp.mapplot()

    # Compute the moments at each position to come up with reasonable guesses.
    # This speeds up the process enormously, but can easily mess up the fits if
    # there are bad pixels
    sp.momenteach(vheight=vheight, verbose=verbose)
    sp.fiteach(errmap=errmap, multifit=None, verbose_level=verbose_level,
               signal_cut=signal_cut, usemomentcube=True, blank_value=np.nan,
               verbose=verbose, **kwargs)

    # steal the header from the error map
    f = fits.open(cubefilename)
    # start replacing components of the fits object
    f[0].data = np.concatenate([sp.parcube,sp.errcube,sp.integralmap])
    f[0].header['PLANE1'] = 'amplitude'
    f[0].header['PLANE2'] = 'velocity'
    f[0].header['PLANE3'] = 'sigma'
    f[0].header['PLANE4'] = 'err_amplitude'
    f[0].header['PLANE5'] = 'err_velocity'
    f[0].header['PLANE6'] = 'err_sigma'
    f[0].header['PLANE7'] = 'integral'
    f[0].header['PLANE8'] = 'integral_error'
    f[0].header['CDELT3'] = 1
    f[0].header['CTYPE3'] = 'FITPAR'
    f[0].header['CRVAL3'] = 0
    f[0].header['CRPIX3'] = 1
    # save your work
    if astropy.version.major >= 2 or (astropy.version.major==1 and astropy.version.minor>=3):
        f.writeto(outfilename, overwrite=overwrite)
    else:
        f.writeto(outfilename, clobber=overwrite)

    return sp

