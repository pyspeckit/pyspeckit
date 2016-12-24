Fitting H2 1-0 S(1) line in NIR data cube
=========================================

.. code-block:: python

  import astropy
  import pyspeckit
  import pylab as pl
  import numpy as np
  import astropy.io.fits as pyfits


  # load the cube as a pyspeckit Cube
  # I usually get an error here if cube.fits' header doesn't have 'CTYPE3' = 'WAVE'
  cube = pyspeckit.Cube('cube.fits')

  # Slice the cube over the wavelength range you'd like to fit
  cube_h2 = cube.slice(21145,21250,unit='Angstrom')

  # rescale data to make fitting nicer
  # I add an offset to avoid negative values
  cube_h2.cube *= 1e17 
  cube_h2.cube += 30.

  # Do an initial plot & fit of a single spectrum
  # at a pixel with good S/N
  cube_h2.plot_spectrum(100,100)

  ### The following command is just for setup!  The actual fitting occurs below. ###
  # Here I'm fitting two gaussians with 4 parameters each (background offset, 
  # amplitude, wavelength centroid, linewidth).
  # I find that if I let both backgrounds be free parameters, pyspeckit returns
  # unrealistic values for both backgrounds, so I fix the 2nd gaussian's background
  # level to 0.  The actual command to fix the parameter comes in the fiteach call.
  cube_h2.specfit(fittype='vheightgaussian',guesses=[3,1,2.1220e4,2,0,50,21232.,2],quiet=False,save=False)

  # Get ready for the interactive plots that come up after fiteach finishes 
  cube_h2.mapplot.makeplane(estimator=np.nansum)

  # For my cube.fits, (21145,21195) covers a part of the spectrum that is free of
  # spectral lines.  The std variable will be used to estimate the S/N of a line 
  # during fiteach.
  std = cube_h2.stats((21145,21195))['std']

  #### Here's where all the fitting happens.
  ## With the "parlimited" and "parlimits" keywords, I have restricted
  ## the range for the wavelength centroid and linewidth parameters.
  ## With the "fixed" keyword, I have held the 2nd gaussian's background level
  ## to zero, and the "signal_cut" keyword rejects fits for voxels below a 
  ## user-specified S/N threshold.
  cube_h2.fiteach(use_nearest_as_guess=False,
                  guesses=[3,1,2.1220e4,2,0,50,21232.,2],
                  fittype='vheightgaussian',
                  integral=False,
                  multicore=4,
                  negamp=False,
                  verbose_level=2,
                  errspec=np.ones(cube_h2.shape)*std,
                  parlimited=[(False,False), (False,False), (True,True),
                              (True,True), (False,False), (True,True),
                              (True,True), (True,True)],
                  parlimits=[(0.9,1.4), (0,16), (21210,21225), (0.5,5),
                             (0.9,1.4), (0,100), (21227.,21236), (0.5,5)],
                  fixed=[False, False, False, False, True, False, False,
                         False], 
                  signal_cut=20,
                  start_from_point=(100,100))

  # plot the fits as images (you can click on background image to see the spectra + fits)
  amp_max = np.max(cube_h2.parcube[1,:,:])
  cube_h2.mapplot(estimator=1,vmax=amp_max,vmin=0)
  cube_h2.mapplot.axis.set_title("Amplitude")

  cube_h2.mapplot.figure=pl.figure(5)
  cube_h2.mapplot(estimator=3, vmax=5, vmin=0)
  cube_h2.mapplot.axis.set_title("Line Width")

  cube_h2.mapplot.figure=pl.figure(6)
  cube_h2.mapplot(estimator=2,vmin=21215,vmax=21225)
  cube_h2.mapplot.axis.set_title("Line Center")

  cube_h2.mapplot.figure=pl.figure(7)
  cube_h2.mapplot(estimator=0,vmax=100,vmin=0)
  cube_h2.mapplot.axis.set_title("Background")
  pl.show()

  ## Create the image
  background = (cube_h2.parcube[0,:,:] - 30.) / 1e17
  amplitude = cube_h2.parcube[1,:,:] / 1e17
  linecenter = cube_h2.parcube[2,:,:]
  sigma = cube_h2.parcube[3,:,:] / 3e5 * h2_linecenter
  image = np.sqrt(2*np.pi)*h2_amplitude*h2_sigma

  # Clean up the header  
  # (this is a bit of a hacky way to do it, but it works)
  cube_h2.header['NAXIS'] = 2
  del cube_h2.header['NAXIS3']
  # a nicer way is to use WCS:
  from astropy import wcs
  newheader = wcs.WCS(cube_h2.header).sub([wcs.WCSSUB_CELESTIAL]).to_header()
  cube2.header = newheader
  # however, this approach may lose other important header keywords


  # Write the image to file
  h2filename = input_filename.replace("cube.fits","h2_1-0S1.fits")
  h2fits = pyfits.PrimaryHDU(data=h2_image,header=cube_h2.header)
  if astropy.version.major >= 2 or (astropy.version.major==1 and astropy.version.minor>=3):
       h2fits.writeto(h2filename,overwrite=True)
  else:
       h2fits.writeto(h2filename,clobber=True)
       
  # Write pyspeckit parcube and errcube to file
  pyspeckit_fits_filename = input_filename.replace("cube.fits",
                                                   "pyspeckitfits_h2_1-0S1.fits")
  cube_h2.write_fit(pyspeckit_fits_filename,overwrite=True)

