import pyspeckit
import pyfits
import numpy as np
import os

def cube_fit(cubefilename, outfilename, errfilename=None, scale_keyword=None,
        vheight=False, verbose=False, signal_cut=3, verbose_level=2,
        clobber=True, **kwargs):
    """
    Light-weight wrapper for cube fitting

    Takes a cube and error map (error will be computed naively if not given)
    and computes moments then fits for each spectrum in the cube.  It then
    saves the fitted parameters to a reasonably descriptive output file.

    **kwargs are passed to pyspeckit.Spectrum.specfit
    """

    # Load the spectrum
    sp = pyspeckit.Cube(cubefilename,scale_keyword=scale_keyword)
    if os.path.exists(errfilename):
        errmap = pyfits.getdata(errfilename)
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
    sp.fiteach(errmap=errmap, multifit=True, verbose_level=verbose_level,
            signal_cut=signal_cut, usemomentcube=True, blank_value=np.nan,
            verbose=verbose, **kwargs)

    # steal the header from the error map
    f = pyfits.open(cubefilename)
    # start replacing components of the pyfits object
    f[0].data = np.concatenate([sp.parcube,sp.errcube,sp.integralmap])
    f[0].header.update('PLANE1','amplitude')
    f[0].header.update('PLANE2','velocity')
    f[0].header.update('PLANE3','sigma')
    f[0].header.update('PLANE4','err_amplitude')
    f[0].header.update('PLANE5','err_velocity')
    f[0].header.update('PLANE6','err_sigma')
    f[0].header.update('PLANE7','integral')
    f[0].header.update('PLANE8','integral_error')
    f[0].header.update('CDELT3',1)
    f[0].header.update('CTYPE3','FITPAR')
    f[0].header.update('CRVAL3',0)
    f[0].header.update('CRPIX3',1)
    # save your work
    f.writeto(outfilename, clobber=clobber)

    return sp

