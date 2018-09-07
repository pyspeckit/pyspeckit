from __future__ import print_function
import pyspeckit

# Grab a .fits spectrum with a legitimate header
sp = pyspeckit.Spectrum('G031.947+00.076_nh3_11_Tastar.fits')
""" HEADER:
SIMPLE  =                    T / Written by IDL:  Tue Aug 31 18:17:01 2010
BITPIX  = -64
NAXIS   =                    1 / number of array dimensions
NAXIS1  =                 8192 /Number of positions along axis 1
CDELT1  =  -0.077230503
CRPIX1  =  4096.0000
CRVAL1  =  68.365635
CTYPE1  =  'VRAD'
CUNIT1  =  'km/s    '
SPECSYS =  'LSRK'
RESTFRQ =  2.3694500e+10
VELOSYS =  -43755.930
CDELT1F =  6103.5156
CRPIX1F =  4096.0000
CRVAL1F =  2.3692555e+10
CTYPE1F =  'FREQ'
CUNIT1F =  'Hz'
SPECSYSF=  'LSRK'
RESTFRQF=  2.3694500e+10
VELOSYSF=  -43755.930
VDEF    =  'RADI-LSR'
SRCVEL  =  70.000000
ZSOURCE  =  0.00023349487
BUNIT   = 'K       '
OBJECT  =  'G031.947+00.076'
TELESCOP=  'GBT'
TSYS    =  42.1655
ELEV    =  34.904846
AIRMASS =  1.7475941
LINE    =  'nh3_11'
FREQ    =  23.692555
TARGLON =  31.947236
TARGLAT =  0.076291610
MJD-AVG =  54548.620
CONTINUU=  0.0477613
CONTERR =  0.226990
SMTHOFF =  0
COMMENT   1 blank line
END
"""

# Start by computing the error using a reasonably signal-free region
sp.error[:] = sp.stats((-100, 50))['std']

# Change the plot range to be a reasonable physical coverage (the default is to
# plot the whole 8192 channel spectrum)
sp.plotter(xmin=-100,xmax=300)

# There are many extra channels, so let's smooth.  Default is a Gaussian
# smooth.  Downsampling helps speed up the fitting (assuming the line is still
# Nyquist sampled, which it is)
sp.smooth(2)
# replot after smoothing
sp.plotter(xmin=-100,xmax=300)

# First, fit a gaussian to the whole spectrum as a "first guess" (good at
# picking up the centroid, bad at getting the width right)
# negamp=False forces the fitter to search for a positive peak, not the
# negatives created in this spectrum by frequency switching
sp.specfit.selectregion(xmin=60,xmax=120,xtype='wcs')
sp.specfit(negamp=False, guesses='moments')
# Save the fit...
sp.plotter.figure.savefig('nh3_gaussfit.png')
# and print some information to screen
print("Guesses: ", sp.specfit.guesses)
print("Best fit: ", sp.specfit.modelpars)

# Run the ammonia spec fitter with a reasonable guess
# Since we only have a single line (1-1), the kinetic temperature is
# unconstrained: we'll fix it at 7 K.  Similarly, the ortho fraction
# is fixed to 0.5
T=True; F=False
sp.specfit(fittype='ammonia_tau',
           guesses=[7,4.45,4.5,0.84,96.2,0.43],
           fixed=[T,F,F,F,F,T],
           quiet=False)

# plot up the residuals in a different window.  The residuals strongly suggest
# the presence of a second velocity component.
sp.specfit.plotresiduals()

sp.plotter.figure.savefig('nh3_ammonia_vtau_fit.png')
print("Guesses: ", sp.specfit.guesses)
print("Best fit: ", sp.specfit.modelpars)


# re-plot zoomed in
sp.plotter(xmin=70,xmax=125)
# replot the fit
sp.specfit.plot_fit()
sp.plotter.figure.savefig('nh3_ammonia_fit_vtau_zoom.png')

# refit with two components
sp.specfit(fittype='ammonia_tau',
           guesses=[7,3.5,4.5,0.68,97.3,0.5]+[7,4.2,4.5,0.52,95.8,0.5],
           fixed=[T,F,F,F,F,T]*2,
           quiet=False)
sp.specfit.plotresiduals()
sp.plotter.figure.savefig('nh3_ammonia_multifit_vtau_zoom.png')

# compare to the 'thin' version
# In the thin version, Tex = Tk by force
sp.specfit.Registry.add_fitter('ammonia_tau_thin',
                               pyspeckit.spectrum.models.ammonia.ammonia_model_vtau_thin(),
                               5)
sp.specfit(fittype='ammonia_tau_thin',
           guesses=[7,4.5,0.68,97.3,0.5]+[7,4.5,0.52,95.8,0.5],
           fixed=[F,F,F,F,T]*2,
           quiet=False)
sp.specfit.plotresiduals()
sp.plotter.figure.savefig('nh3_ammonia_multifit_vtau_thin_zoom.png')
