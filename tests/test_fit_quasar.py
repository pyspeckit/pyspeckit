from matplotlib import pylab
pylab.ion()
C_IV=1549.48
He_II=1640.4
Si_IV=1397.61
Si_O_IV=1399.8
guesses = [0.1, C_IV, 3, 
           0.05, C_IV, 10,
           0.02, He_II, 20,
           0.08, Si_IV, 3,
           0.08, Si_O_IV, 3]
import pyspeckit
sp = pyspeckit.Spectrum('spectrum_4.txt')
sp.downsample(10)
sp.data *= 1e15
sp.error[:] = 1 # just a guess...

# set up units properly
sp.xarr.units='angstroms'
sp.xarr.xtype = 'wavelength'
sp.units = r'$1\times10^{-15}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$'


# plot!
sp.plotter(xmin=1350,xmax=1820)
sp.plotter(xmin=1000,xmax=3520,ymax=1)

sp.baseline.powerlaw=True
sp.baseline(xmin=1200, xmax=3500,exclude=[1200,1440,1490,1700,1830,1995,2700,2850],subtract=False,
        reset_selection=True, highlight_fitregion=True,
        powerlaw=True,quiet=False,LoudDebug=True)
# uncomment all of this if you want to re-scale back down by 1e15
#sp.baseline.baselinepars[0] /= 1e15
#sp.data /= 1e15
#sp.error /= 1e15
#sp.baseline.set_basespec_frompars()
#sp.baseline.set_spectofit()
#sp.baseline.plot_baseline()
#sp.baseline.highlight_fitregion()
pylab.show()

# crop for faster fitting
small_sp = sp.copy()
small_sp.crop(1000,3520)
small_sp.plotter(xmin=1000,xmax=3520,ymax=1)
small_sp.baseline.powerlaw=True
small_sp.baseline(xmin=1200, xmax=3500,exclude=[1200,1440,1490,1700,1830,1995,2700,2850],subtract=False,
        reset_selection=True, highlight_fitregion=True,
        powerlaw=True,quiet=False,LoudDebug=True)

small_sp.plotter(xmin=1300,xmax=1700,ymax=1)
small_sp.specfit(guesses=guesses, annotate=False, fittype='gaussian', quiet=False)
