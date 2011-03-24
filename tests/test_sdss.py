import spectrum
spec = spectrum.Spectrum('sample_sdss.txt',errspecnum=2)
spec.plotter(xmin = 6450, xmax = 6675, ymin = -20, ymax = 120)
print "Plotter min/max: ",spec.plotter.xmin,spec.plotter.xmax," Fitter min/max: ",spec.specfit.gx1,spec.specfit.gx2," Fitregion= ",spec.baseline.fitregion
spec.specfit(guesses = [55, 30, 6800, 60], annotate = False, interactive = False, quiet=False)
print "Plotter min/max: ",spec.plotter.xmin,spec.plotter.xmax," Fitter min/max: ",spec.specfit.gx1,spec.specfit.gx2," Fitregion= ",spec.baseline.fitregion
raw_input("Wait here a moment")
spec.specfit(guesses = [50, 6549.86, 5, 100, 6564.614, 5, 50, 6564.614, 50, 50, 6585.27, 5], tied = ['', '', '', '', '', '', '', '', '', '3 * p[0]', '', 'p[2]'])
raw_input("Wait here a moment")
spec.specfit.plotresiduals()
raw_input("Wait here a moment")
