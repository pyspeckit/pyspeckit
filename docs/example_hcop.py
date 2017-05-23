import pyspeckit

# load a FITS-compliant spectrum
spec = pyspeckit.Spectrum('10074-190_HCOp.fits')
# The units are originally frequency (check this by printing spec.xarr.units).
# I want to know the velocity.  Convert!
# Note that this only works because the reference frequency is set in the header
spec.xarr.convert_to_unit('km/s')

# plot it
spec.plotter()

# Subtract a baseline
spec.baseline()

# Fit a gaussian.  We know it will be an emission line, so we force a positive guess
spec.specfit(negamp=False)

# Note that the errors on the fits are larger than the fitted parameters.
# That's because this spectrum did not have an error assigned to it. 
# Let's use the residuals:
spec.specfit.plotresiduals()

# Now, refit with error determined from the residuals:
# (we pass in guesses to save time / make sure nothing changes)
spec.specfit(guesses=spec.specfit.modelpars)

# Save the figures to put on the web....
spec.plotter.figure.savefig("simple_fit_example_HCOp.png")
spec.specfit.residualaxis.figure.savefig("simple_fit_example_HCOp_residuals.png")

# Also, let's crop out stuff we don't want...
spec.crop(-100,100)
# replot after cropping (crop doesn't auto-refresh)
spec.plotter()
# replot the fit without re-fitting
spec.specfit.plot_fit()
# show the annotations again
spec.specfit.annotate()

spec.plotter.figure.savefig("simple_fit_example_HCOp_cropped.png")
