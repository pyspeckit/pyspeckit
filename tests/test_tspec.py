from pyspeckit import spectrum
tsp = spectrum.Spectrum('BD60596_merged.fits',filetype='tspec')
tsp.plotter()
tsp.baseline(include=[0.96,0.99,1.1655735969543457, 1.189874529838562,
    1.206290602684021, 1.263965368270874, 1.3026881217956543,
    1.3405940532684326,1.5,1.78,2.03,2.46],powerlaw=True,subtract=False,quiet=False,LoudDebug=True)
tsp.baseline.show_fitregion()
#tsp.plotter(xmin=1.20,xmax=1.35,reset_ylimits=True,clear=False)
#tsp.plotter(xmin=1.260,xmax=1.303,errstyle='fill',figure=1)
#tsp.baseline(interactive=True,powerlaw=True,quiet=False,subtract=False)
#tsp.baseline(exclude=[1.2779,1.2843],order=2)

#tsp.specfit(guesses=[-1.95e-14, 1.2800, 1e-3, 1e-4, 0.5e-14, 1.2801,1e-3, 1e-4],quiet=False,shh=False,fittype='voigt')

"""
This is all old debug stuff
#tsp.plotter(reset_ylimits=True)
#tsp.specfit(interactive=True)
#tsp.baseline(interactive=True)

tsp2 = spectrum.Spectrum('BD60596_merged.fits',filetype='tspec')
tsp2.crop(1.260,1.305)
tsp2.plotter(errstyle='fill')
tsp2.baseline(exclude=[1.2779,1.2843],order=2)
tsp2.specfit(guesses=[-1.95e-14, 1.2800, 1e-3, 0.5e-14, 1.2801,1e-3],quiet=False,shh=False)

tsp3 = spectrum.Spectrum('BD60596_merged.fits',filetype='tspec')
tsp3.crop(1.260,1.305)
tsp3.data *= 1e14
tsp3.error *= 1e14
tsp3.plotter(errstyle='fill')
tsp3.baseline(exclude=[1.2779,1.2843],order=2)
tsp3.specfit(guesses=[-1.95, 1.2800, 1e-3, 0.5, 1.2801,1e-3],quiet=False,shh=False)

tsp4 = spectrum.Spectrum('BD60596_merged.fits',filetype='tspec')
tsp4.crop(1.260,1.305)
tsp4.xarr.convert_to_unit('angstroms')
tsp4.data *= 1e14
tsp4.error *= 1e14
tsp4.plotter(errstyle='fill')
tsp4.baseline(exclude=[1.2779e4,1.2843e4],order=4)
tsp4.specfit(guesses=[-1.7879420985498577, 1.2808937741295388e4, 0.0012223729309375341e4, 3.0317627598183574e-1, 1.2813118500364311e4, 0.0011904119556266335e4],quiet=False,shh=False)
tsp4.specfit(guesses=[-0.95, 1.2826e4, 10, -0.8, 1.2792e4,5],quiet=False,shh=False)

import pylab
pylab.show()
"""
import pylab as pl
pl.show()
