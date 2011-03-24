import spectrum
spec = spectrum.Spectrum('sample_sdss.fits',errspecnum=2)
spec.plotter(xmin = 6500, xmax = 7200, ymin = 40, ymax = 120)
print "Plotter min/max: ",spec.plotter.xmin,spec.plotter.xmax," Fitter min/max: ",spec.specfit.gx1,spec.specfit.gx2," Fitregion= ",spec.baseline.fitregion
spec.specfit(guesses = [55, 30, 6800, 60], annotate = False, interactive = False, quiet=False)
print "Plotter min/max: ",spec.plotter.xmin,spec.plotter.xmax," Fitter min/max: ",spec.specfit.gx1,spec.specfit.gx2," Fitregion= ",spec.baseline.fitregion
raw_input("Wait here a moment")
spec.specfit(guesses = [30, 6800, 60, 5, 6863, 3, 5, 6883, 3], annotate = False, interactive = False, quiet=False)
raw_input("Wait here a moment")
spec.specfit(guesses = [30, 6800, 60, 5, 6863, 3, 5, 6883, 3, 5, 7010, 3, 5, 7040, 3], annotate = False, interactive = False, quiet=False)
raw_input("Wait here a moment")
spec.specfit.plotresiduals()


raw_input("Wait here a moment")
