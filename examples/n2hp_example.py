import pyspeckit

sp = pyspeckit.Spectrum('n2hp_opha_example.fits')

sp.Registry.add_fitter('n2hp_vtau',pyspeckit.models.n2hp.n2hp_vtau_fitter,4,multisingle='multi')

sp.specfit(fittype='n2hp_vtau',multifit=True,guesses=[15,2,4,0.2])

sp.plotter()
sp.specfit(fittype='n2hp_vtau',multifit=True,guesses=[15,2,4,0.2],show_components=True)
