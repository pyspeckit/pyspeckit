from pyspeckit import spectrum

if not 'interactive' in globals():
    interactive=False

spec = spectrum.Spectrum('sample_sdss.txt',errspecnum=2)
spec.plotter(xmin = 6450, xmax = 6775, ymin = -20, ymax = 120)
print "Plotter min/max: ",spec.plotter.xmin,spec.plotter.xmax," Fitter min/max: ",spec.specfit.gx1,spec.specfit.gx2," Fitregion= ",spec.baseline.fitregion
if interactive: raw_input("Wait here a moment") 
NIIa = 6549.86
NIIb = 6585.27
Halpha = 6564.614
SIIa = 6718.29
SIIb = 6732.68
# Use [SII] lines to model narrow lines, then force [NII] lines and narrow H-alpha to have same width as [SII].  
# Fit 2 additional broad components to H-alpha.
spec.specfit(guesses = [50, NIIa, 5, 100, Halpha, 5, 50, Halpha, 50, 50, NIIb, 5, 20, SIIa, 5, 20, SIIb, 5, 20, Halpha, 5], 
    tied = ['', '', 'p[17]', '', '', 'p[17]', '', 'p[4]', '', '3 * p[0]', '', 'p[17]', '', '', 'p[17]', '','','', '', '', ''])
#if interactive: raw_input("Wait here a moment")
#spec.specfit.plotresiduals()
#if interactive: raw_input("Wait here a moment")
raw_input('Done.')
