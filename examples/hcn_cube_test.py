import pyspeckit
try:
    from astropy.io import fits as pyfits
except ImportError:
    import pyfits
import numpy as np

# Load the spectrum
sp = pyspeckit.Cube('region5_hcn_crop.fits')
errmap = pyfits.getdata('region5.hcn.errmap.fits')

# Register the fitter
# The N2H+ fitter is 'built-in' but is not registered by default; this example
# shows how to register a fitting procedure
# 'multi' indicates that it is possible to fit multiple components and a background will not automatically be fit
# 4 is the number of parameters in the model (excitation temperature, optical depth, line center, and line width)
sp.Registry.add_fitter('hcn_amp',pyspeckit.models.hcn.hcn_amp,3)

# Run the fitter
sp.mapplot()

# use an individual spectrum selected semi-arbitrarily from the map to get an estimate of the error
# this method has been rendered obsolete - use the error map instead
#s = sp.get_spectrum(20,20)
#s.plotter()
#s.Registry = sp.Registry
#s.specfit.Registry = sp.Registry
#s.specfit(fittype='hcn_amp',guesses=[2.5,-5.6,1.5],show_components=True,debug=True,quiet=False)
#s.specfit(fittype='hcn_amp',guesses=[2.5,-5.6,1.5],show_components=True,debug=True,quiet=False)
#sp.error = s.specfit.errspec


# Compute the moments at each position to come up with reasonable guesses.
# This speeds up the process enormously, but can easily mess up the fits if
# there are bad pixels
sp.momenteach(vheight=False, verbose=False)
sp.momentcube[2,:,:] /= 2.5 # the HCN line profile makes the fitter assume a 2.5x too large line
sp.fiteach(fittype='hcn_amp', errmap=errmap,
        guesses=[1.0,-5.6,1.5], verbose_level=2, signal_cut=4,
        usemomentcube=True, blank_value=np.nan, verbose=False,
        direct=True, multicore=4)

# steal the header from the error map
f = pyfits.open('region5.hcn.errmap.fits')
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
f.writeto('region5.hcn.nosmooth.fit.fits', clobber=True)

