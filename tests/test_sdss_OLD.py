import pyspeckit

if not 'interactive' in globals():
    interactive=False
if not 'savedir' in globals():
    savedir = ''

# Rest wavelengths of the lines we are fitting - use as initial guesses
NIIa = 6549.86
NIIb = 6585.27
Halpha = 6564.614
SIIa = 6718.29
SIIb = 6732.68

# Initialize spectrum object and plot region surrounding Halpha-[NII] complex
spec = pyspeckit.Spectrum('sample_sdss.txt', errorcol=2)
#spec.units = 'erg s^{-1} cm^{-2} \\AA^{-1}'
#spec.xarr.units='angstroms'
spec.xarr.redshift=0.05
spec.crop(6450,6755)
spec.plotter(xmin = 6450, xmax = 6775, ymin = 0, ymax = 150)

# Use [SII] lines to model narrow lines, then force [NII] lines and narrow H-alpha to have same width as [SII].  
# Will fit 1 additional broad component to H-alpha (standard for AGN spectra)
# Wavelengths are all tied together
guesses = [50, NIIa, 5, 100, Halpha, 5, 50, Halpha, 50, 50, NIIb, 5, 20, SIIa, 5, 20, SIIb, 5]
tied = ['', '', 'p[17]', '', '', 'p[17]', '', 'p[4]', '', '3 * p[0]', '', 'p[17]', '', '', 'p[17]', '', '', '']

# Actually do the fit.
spec.specfit(multifit=True, guesses=guesses, tied=tied, annotate=False, show_components=True, debug=True, quiet=False)
spec.plotter.refresh()

# Let's use the measurements class to derive information about the emission lines.
spec.measure(z = 0.05, fluxnorm = 1e-17, debug=True)

# Now overplot positions of lines and annotate
y = spec.plotter.ymax * 0.85
for i, line in enumerate(spec.measurements.lines.keys()):
    
    if line not in spec.speclines.optical.lines.keys(): continue
    
    x = spec.measurements.lines[line]['modelpars'][1]
    spec.plotter.axis.plot([x]*2, [spec.plotter.ymin, spec.plotter.ymax], ls = '--', color = 'k')
    
    spec.plotter.axis.annotate(spec.speclines.optical.lines[line][-1], (x, y), rotation = 90, ha = 'right', va = 'center')

# Make some nice axis labels
spec.plotter.axis.set_xlabel(r'Wavelength $(\AA)$')
spec.plotter.axis.set_ylabel(r'Flux $(10^{-17} \mathrm{erg/s/cm^2/\AA})$')
spec.plotter.refresh()

# Print out spectral line information
print "Line   Flux (erg/s/cm^2)     Amplitude (erg/s/cm^2)    FWHM (Angstrom)   Luminosity (erg/s)"
for line in spec.measurements.lines.keys():
    print line, spec.measurements.lines[line]['flux'], spec.measurements.lines[line]['amp'], \
        spec.measurements.lines[line]['fwhm'], spec.measurements.lines[line]['lum']

spec.specfit.add_sliders()
        
# Notice that because we supplied the objects redshift and flux normalization, the measurements class
# automatically calculated line luminosities.  Also, it separates the broad and narrow H-alpha components, and identifies which lines are which. How nice!

#raw_input("Done.")

# unnecessary?
#spec.specfit.plot_fit(show_components=True)

#raw_input("Done (again).")

# Save the figure
spec.plotter.savefig(savedir+"sdss_fit_example.png")

# Do the same thing using lineid_plot
try: 
    import lineid_plot
    import pylab

    spec.plotter(figure=pylab.figure())
    spec.specfit.plot_fit(annotate=False)
    spec.plotter.line_ids_from_measurements()

    spec.plotter.savefig(savedir+"sdss_fit_example_lineidplot.png")

except ImportError:
    pass
