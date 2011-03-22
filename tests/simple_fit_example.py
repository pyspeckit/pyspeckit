import spectrum

# load a FITS-compliant spectrum
spec = spectrum.Spectrum('10074-190_HCOp.fits')
# The units are originally frequency (check this by printing spec.xarr.units).
# I want to know the velocity.  Convert!
# Note that this only works because the reference frequency is set in the header
spec.xarr.frequency_to_velocity()
# Default conversion is to m/s, but we traditionally work in km/s
spec.xarr.convert_to_unit('km/s')
# plot it up!
spec.plotter()
# Subtract a baseline (the data is only 'mostly' reduced)
spec.baseline()
# Fit a gaussian.  We know it will be an emission line, so we force a positive guess
spec.specfit(negamp=False)
# Save the figure to put on the web....
spec.plotter.figure.savefig("simple_fit_example_HCOp.png")
