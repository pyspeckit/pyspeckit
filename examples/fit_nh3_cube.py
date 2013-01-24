"""
Fit NH3 Cube
============

Example script to fit all pixels in an NH3 data cube
"""
import pyspeckit
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
import numpy as np
import os
import agpy

# set up CASA-like shortcuts
F=False; T=True

# Some optional parameters for the script
fitcube = True

# Mask out low S/N pixels (to speed things up)
mask = pyfits.getdata('hotclump_11_mask.fits')
mask = np.isfinite(mask) * (mask > 0)

# Load the data using a mask
# Then calibrate the data (the data we're loading in this case are in Janskys,
# but we want surface brightness in Kelvin for the fitting process)
cube11 = pyspeckit.Cube('hotclump_11.cube_r0.5_rerun.image.fits', maskmap=mask)
cube11.cube *= 13.6 * (300.0 / (pyspeckit.spectrum.models.ammonia.freq_dict['oneone']/1e9))**2 * 1./cube11.header.get('BMAJ')/3600. * 1./cube11.header.get('BMIN')/3600. 
cube11.units = "K"
cube22 = pyspeckit.Cube('hotclump_22.cube_r0.5_contsub.image.fits', maskmap=mask)
cube22.cube *= 13.6 * (300.0 / (pyspeckit.spectrum.models.ammonia.freq_dict['twotwo']/1e9))**2 * 1./cube22.header.get('BMAJ')/3600. * 1./cube22.header.get('BMIN')/3600. 
cube22.units = "K"
cube44 = pyspeckit.Cube('hotclump_44.cube_r0.5_contsub.image.fits', maskmap=mask)
cube44.cube *= 13.6 * (300.0 / (pyspeckit.spectrum.models.ammonia.freq_dict['fourfour']/1e9))**2 * 1./cube44.header.get('BMAJ')/3600. * 1./cube44.header.get('BMIN')/3600. 
cube44.units = "K"

# Compute an error map.  We use the 1-1 errors for all 3 because they're
# essentially the same, but you could use a different error map for each
# frequency
errmap11 = pyfits.getdata('hotclump_11.cube_r0.5_rerun.image.moment_linefree.fits').squeeze() * 13.6 * (300.0 / (pyspeckit.spectrum.models.ammonia.freq_dict['oneone']/1e9))**2 * 1./cube11.header.get('BMAJ')/3600. * 1./cube11.header.get('BMIN')/3600. 
errmap11[errmap11 != errmap11] = agpy.smooth(errmap11, interp_nan=True)[errmap11 != errmap11]

# Stack the cubes into one big cube.  The X-axis is no longer linear: there
# will be jumps from 1-1 to 2-2 to 4-4.  
cubes = pyspeckit.CubeStack([cube11,cube22,cube44], maskmap=mask)
cubes.units = "K"

# Make a "moment map" to contain the initial guesses
# If you've already fit the cube, just re-load the saved version
# otherwise, re-fit it
if os.path.exists('hot_momentcube.fits'):
    momentcubefile = pyfits.open('hot_momentcube.fits')
    momentcube = momentcubefile[0].data
else:
    cube11.mapplot()
    # compute the moment at each pixel
    cube11.momenteach()
    momentcube = cube11.momentcube
    momentcubefile = pyfits.PrimaryHDU(data=momentcube, header=cube11.header)
    momentcubefile.writeto('hot_momentcube.fits',clobber=True)

# Create a "guess cube".  Because we're fitting physical parameters in this
# case, we want to make the initial guesses somewhat reasonable
# As above, we'll just reload the saved version if it exists
guessfn = 'hot_guesscube.fits'
if os.path.exists(guessfn):
    guesscube = pyfits.open(guessfn)
    guesses = guesscube[0].data
else:
    guesses = np.zeros((6,)+cubes.cube.shape[1:])
    guesses[0,:,:] = 20                    # Kinetic temperature 
    guesses[1,:,:] = 5                     # Excitation  Temp
    guesses[2,:,:] = 14.5                  # log(column)
    guesses[3,:,:] = momentcube[3,:,:] / 5 # Line width / 5 (the NH3 moment overestimates linewidth)                  
    guesses[4,:,:] = momentcube[2,:,:]     # Line centroid              
    guesses[5,:,:] = 0.5                   # F(ortho) - ortho NH3 fraction (fixed)

    guesscube = pyfits.PrimaryHDU(data=guesses, header=cube11.header)
    guesscube.writeto(guessfn, clobber=True)

# This bit doesn't need to be in an if statment
if fitcube:
    # excise guesses that fall out of the "good" range
    guesses[4,:,:][guesses[4,:,:] > 100] = 100.0
    guesses[4,:,:][guesses[4,:,:] < 91] = 95

    # do the fits
    # signal_cut means ignore any pixel with peak S/N less than this number
    # In this fit, many of the parameters are limited
    # start_from_point selects the pixel coordinates to start from
    # use_nearest_as_guess says that, at each pixel, the input guesses will be
    # set by the fitted parameters from the nearest pixel with a good fit
    # HOWEVER, because this fitting is done in parallel (multicore=12 means
    # 12 parallel fitting processes will run), this actually means that EACH
    # core will have its own sub-set of the cube that it will search for good 
    # fits. So if you REALLY want consistency, you need to do the fit in serial.
    cubes.fiteach(fittype='ammonia', multifit=True, guesses=guesses,
            integral=False, verbose_level=3, fixed=[F,F,F,F,F,T], signal_cut=3,
            limitedmax=[F,F,F,F,T,T],
            maxpars=[0,0,0,0,101,1],
            limitedmin=[T,T,F,F,T,T],
            minpars=[2.73,2.73,0,0,91,0],
            use_nearest_as_guess=True, start_from_point=(94,250),
            multicore=12,
            errmap=errmap11)

    # Save the fitted parameters in a data cube
    fitcubefile = pyfits.PrimaryHDU(data=np.concatenate([cubes.parcube,cubes.errcube]), header=cubes.header)
    fitcubefile.header.update('PLANE1','TKIN')
    fitcubefile.header.update('PLANE2','TEX')
    fitcubefile.header.update('PLANE3','COLUMN')
    fitcubefile.header.update('PLANE4','SIGMA')
    fitcubefile.header.update('PLANE5','VELOCITY')
    fitcubefile.header.update('PLANE6','FORTHO')
    fitcubefile.header.update('PLANE7','eTKIN')
    fitcubefile.header.update('PLANE8','eTEX')
    fitcubefile.header.update('PLANE9','eCOLUMN')
    fitcubefile.header.update('PLANE10','eSIGMA')
    fitcubefile.header.update('PLANE11','eVELOCITY')
    fitcubefile.header.update('PLANE12','eFORTHO')
    fitcubefile.header.update('CDELT3',1)
    fitcubefile.header.update('CTYPE3','FITPAR')
    fitcubefile.header.update('CRVAL3',0)
    fitcubefile.header.update('CRPIX3',1)
    fitcubefile.writeto("hot_fitcube_try6.fits")
else: # you can read in a fit you've already done!
    cubes.load_model_fit('hot_fitcube_try6.fits', 6, 'ammonia', _temp_fit_loc=(94,250))
    cubes.specfit.parinfo[5]['fixed'] = True


# Now do some plotting things
import pylab as pl

# Set the map-to-plot to be the line centroid
cubes.mapplot.plane = cubes.parcube[4,:,:]
cubes.mapplot(estimator=None,vmin=91,vmax=101)

# Set the reference frequency to be the 1-1 line frequency
cubes.xarr.refX = pyspeckit.spectrum.models.ammonia.freq_dict['oneone']
cubes.xarr.refX_units='Hz'

# If you wanted to view the spectra in velocity units, use this:
#cubes.xarr.convert_to_unit('km/s')
#cubes.plotter.xmin=55
#cubes.plotter.xmax=135

# Now replace the cube's plotter with a "special" plotter
# The "special" plotter puts the 1-1, 2-2, and 4-4 lines in their own separate
# windows

cubes.plot_special = pyspeckit.wrappers.fitnh3.plotter_override
cubes.plot_special_kwargs = {'fignum':3, 'vrange':[55,135]}
cubes.plot_spectrum(160,99)

# make interactive
pl.ion()
pl.show()

# At this point, you can click on any pixel in the image and see the spectrum
# with the best-fit ammonia profile overlaid.
