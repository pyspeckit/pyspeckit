import pyspeckit

# Load the spectrum
sp = pyspeckit.Spectrum('n2hp_opha_example.fits')

# Register the fitter
# The N2H+ fitter is 'built-in' but is not registered by default; this example
# shows how to register a fitting procedure
# 'multi' indicates that it is possible to fit multiple components and a
# background will not automatically be fit 4 is the number of parameters in the
# model (excitation temperature, optical depth, line center, and line width)
sp.Registry.add_fitter('n2hp_vtau',pyspeckit.models.n2hp.n2hp_vtau_fitter,4)
sp.xarr.velocity_convention = 'radio'
# Run the fitter
sp.specfit(fittype='n2hp_vtau',multifit=None,guesses=[15,2,4,0.2])

# Plot the results
sp.plotter()
# Re-run the fitter (to get proper error bars) and show the individual fit components
sp.specfit(fittype='n2hp_vtau', multifit=None, guesses=[15, 2, 4, 0.2],
           show_hyperfine_components=True)

# Save the figure (this step is just so that an image can be included on the web page)
sp.plotter.savefig('n2hp_ophA_fit.png')
