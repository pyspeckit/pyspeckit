MUSE cube fitting
~~~~~~~~~~~~~~~~~

An example analysis of a `MUSE
<https://www.eso.org/sci/facilities/develop/instruments/muse.html>`_ data cube
using both `spectral_cube <spectral-cube.rtfd.org>`_ and pyspeckit

.. code-block:: python

   from astropy.io import fits
   from astropy import units as u
   import numpy as np
   import spectral_cube
   import pyspeckit
   from spectral_cube.spectral_axis import vac_to_air

   lines = {'OI6300': 6302.046,
            'SIII6313': 6313.8,
            'OI6363': 6365.536,
            'NII6548': 6549.85,
            'HAlpha': 6564.61,
            'HBeta': 4862.69,
            'NII6584': 6585.28,
            'HeI6678': 6679.99556,
            'SII6716': 6718.29,
            'SII6731': 6732.67,
            'pa20': 8394.71,
            'pa19': 8415.63,
            'pa18': 8440.27,
            'pa17': 8469.59,
            'pa16': 8504.83,
            'pa15': 8547.73,
            'pa14': 8600.75,
            'pa13': 8667.40,
            'pa12': 8752.86,
            'pa11': 8865.32,
            'pa10': 9017.8 ,
            'pa9' : 9232.2 ,
            'HeI7067': 7067.138,
            'HeI7283': 7283.355,
            'ArIII7135': 7137.8,
            'ArIII7751': 7753.2,
            'SIII9071':9071.1,
            'NeII5756':5756.24,
            'HeI5877':5877.3,
            'OIII5008': 5008.24,
            'OII4960': 4960.3,
   }

   # use spectral_cube to read the data
   cube = spectral_cube.SpectralCube.read('CUBEec_nall.fits', hdu=1)
   cont = cube.spectral_slab(6380*u.AA, 6500*u.AA).apply_numpy_function(np.mean, axis=0)

   # "slabs" are velocity-cutouts of the cube
   # Cube velocity conversion should use vacuum wavelengths
   slabs = {line:
            cube.with_spectral_unit(u.km/u.s, velocity_convention='optical',
                                    rest_value=wl*u.AA)
                .spectral_slab(-200*u.km/u.s, 250*u.km/u.s)
            for line,wl in lines.items()}

   # Compute 1st moments (intensity-weighted velocity)
   for line,slab in slabs.items():
       print "kms",line
       mom1 = slab.moment1(axis=0)
       mom1.write('moments/moment1_125_{0}_kms.fits'.format(line), overwrite=True)
   mean_moment = np.mean([fits.getdata('moments/moment1_125_{0}_kms.fits'.format(line)) for
                          line in slabs], axis=0)
   hdr = fits.getheader('moments/moment1_125_{0}_kms.fits'.format(line))
   outfilename = 'moments/moment1_125_mean_kms.fits'
   if astropy.version.major >= 2 or (astropy.version.major==1 and astropy.version.minor>=3):
       fits.PrimaryHDU(data=mean_moment, header=hdr).writeto(outfilename, overwrite=True)
   else:
       fits.PrimaryHDU(data=mean_moment, header=hdr).writeto(outfilename, clobber=True)

   # Create a cube of many different lines all in velocity
   # This will allow us to measure a velocity more accurate than is possible
   # with one line alone (assuming all lines have the same velocity) because
   # MUSE undersamples its line-spread-function a little
   newcube_shape = (sum(s.shape[0] for s in slabs.values()),) + slabs.values()[0].shape[1:]
   newcube_spaxis = np.concatenate([s.spectral_axis
                                    for s in slabs.values()]).value*u.km/u.s
   sortvect = newcube_spaxis.argsort()
   sortspaxis = newcube_spaxis[sortvect]

   newcube = np.empty(newcube_shape)

   # normalize each spectrum
   ind = 0
   for ii,slab in enumerate(slabs.values()):
       data = (slab.filled_data[:] - cont) / (slab.sum(axis=0) - cont*slab.shape[0])
       newcube[ind:ind+data.shape[0], :, :] = data
       ind += data.shape[0]

   supercube = newcube[sortvect, :, :]

   # Create a pyspeckit cube so we can then fit a gaussian to each spectrum
   pxarr = pyspeckit.units.SpectroscopicAxis(sortspaxis.value, units='km/s')
   pcube = pyspeckit.Cube(cube=supercube, xarr=pxarr)

   # more cores = more faster
   pcube.fiteach(fittype='gaussian', guesses=[1/np.sqrt(np.pi), 10, 50.0],
                 errmap=np.ones(supercube.shape[1:])/10., multicore=40)

   pcube.write_fit('velocity_fits_125.fits', overwrite=True)
