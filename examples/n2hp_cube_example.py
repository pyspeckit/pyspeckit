import astropy
import pyspeckit
import os
import astropy.units as u
import warnings
from astropy import wcs

if not os.path.exists('n2hp_cube.fit'):
    import astropy.utils.data as aud
    from astropy.io import fits

    try:
        f = aud.download_file('ftp://cdsarc.u-strasbg.fr/pub/cats/J/A%2BA/472/519/fits/opha_n2h.fit')
    except Exception as ex:
        # this might be any number of different timeout errors (urllib2.URLError, socket.timeout, etc)
        # travis-ci can't handle ftp:
        # https://blog.travis-ci.com/2018-07-23-the-tale-of-ftp-at-travis-ci
        print("Failed to download from ftp.  Exception was: {0}".format(ex))
        f = aud.download_file('http://cdsarc.u-strasbg.fr/ftp/cats/J/A+A/472/519/fits/opha_n2h.fit')

    with fits.open(f) as ff:
        ff[0].header['CUNIT3'] = 'm/s'
        for kw in ['CTYPE4','CRVAL4','CDELT4','CRPIX4','CROTA4']:
            if kw in ff[0].header:
                del ff[0].header[kw]
        ff.writeto('n2hp_cube.fit')

# Load the spectral cube cropped in the middle for efficiency
with warnings.catch_warnings():
    warnings.filterwarnings('ignore', category=wcs.FITSFixedWarning)
    spc = pyspeckit.Cube('n2hp_cube.fit')[:,25:29,12:16]
# Set the velocity convention: in the future, this may be read directly from
# the file, but for now it cannot be.
spc.xarr.refX = 93176265000.0*u.Hz
spc.xarr.velocity_convention = 'radio'
spc.xarr.convert_to_unit('km/s')

# Register the fitter
# The N2H+ fitter is 'built-in' but is not registered by default; this example
# shows how to register a fitting procedure
# 'multi' indicates that it is possible to fit multiple components and a
# background will not automatically be fit 4 is the number of parameters in the
# model (excitation temperature, optical depth, line center, and line width)
spc.Registry.add_fitter('n2hp_vtau', pyspeckit.models.n2hp.n2hp_vtau_fitter, 4)

# Get a measurement of the error per pixel
errmap = spc.slice(20, 28, unit='km/s').cube.std(axis=0)

# A good way to write a cube fitter is to have it load from disk if the cube
# fit was completed successfully in the past
if os.path.exists('n2hp_fitted_parameters.fits'):
    spc.load_model_fit('n2hp_fitted_parameters.fits', npars=4, npeaks=1)
else:
    # Run the fitter
    # Estimated time to completion ~ 2 minutes
    spc.fiteach(fittype='n2hp_vtau',
                guesses=[5,0.5,3,1], # Tex=5K, tau=0.5, v_center=12, width=1 km/s
                signal_cut=3, # minimize the # of pixels fit for the example
                start_from_point=(2,2), # start at a pixel with signal
                errmap=errmap,
                use_neighbor_as_guess=True,
                multicore=4,
                )
    # There are a huge number of parameters for the fiteach procedure.  See:
    # http://pyspeckit.readthedocs.org/en/latest/example_nh3_cube.html
    # http://pyspeckit.readthedocs.org/en/latest/cubes.html?highlight=fiteach#pyspeckit.cubes.SpectralCube.Cube.fiteach
    #
# Unfortunately, a complete tutorial on this stuff is on the to-do list;
# right now the use of many of these parameters is at a research level.
# However, pyspeckit@gmail.com will support them!  They are being used
# in current and pending publications

# Save the fitted parameters to a FITS file, and overwrite one if one exists
spc.write_fit('n2hp_fitted_parameters.fits', overwrite=True)

# Show an image of the fitted tex
spc.mapplot(estimator=1, vmin=0, vmax=10)
# you can click on any pixel to see its spectrum & fit

# plot one of the fitted spectra
spc.plot_spectrum(2, 2, plot_fit=True)
# spc.parcube[:,27,14] = [ 14.82569198,   1.77055642,   3.15740051,   0.16035407]
# Note that the optical depth is the "total" optical depth, which is
# distributed among 15 hyperfine components.  You can see this in
# pyspeckit.spectrum.models.n2hp.line_strength_dict
# As a sanity check, you can see that the brightest line has 0.259 of the total
# optical depth, so the peak line brightness is:
# (14.825-2.73) * (1-np.exp(-1.77 * 0.259)) = 4.45
# which matches the peak of 4.67 pretty well

# Show an image of the best-fit velocity
spc.mapplot.plane = spc.parcube[2,:,:]
spc.mapplot(estimator=None)

# running in script mode, the figures won't show by default on some systems
# import pylab as pl
# pl.draw()
# pl.show()
