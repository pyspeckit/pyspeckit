import pyspeckit

# Load up the n2hp model
import pyspeckit.wrappers.n2hp_wrapper
# replace this with the path to your RADEX model
pyspeckit.wrappers.n2hp_wrapper.path_to_radex = '/Users/adam/work/n2hp/'

# Load the spectral cube
spc = pyspeckit.Cube('n2hp_cube.fits')

# Register the fitter
# The N2H+ fitter is 'built-in' but is not registered by default; this example
# shows how to register a fitting procedure
# 'multi' indicates that it is possible to fit multiple components and a
# background will not automatically be fit 4 is the number of parameters in the
# model (excitation temperature, optical depth, line center, and line width)
spc.Registry.add_fitter('n2hp_vtau',pyspeckit.models.n2hp.n2hp_vtau_fitter,4,multisingle='multi')

# Run the fitter
spc.fiteach(fittype='n2hp_vtau', multifit=True,
        guesses=[5,0.5,3,1], # Tex=5K, tau=0.5, v_center=12, width=1 km/s
        )
# There are a huge number of parameters for the fiteach procedure.  See:
# http://pyspeckit.readthedocs.org/en/latest/example_nh3_cube.html
# http://pyspeckit.readthedocs.org/en/latest/cubes.html?highlight=fiteach#pyspeckit.cubes.SpectralCube.Cube.fiteach
#
# Unfortunately, a complete tutorial on this stuff is on the to-do list;
# right now the use of many of these parameters is at a research level.
# However, pyspeckit@gmail.com will support them!  They are being used
# in current and pending publications

# Show an integrated image
spc.mapplot()

# plot one of the fitted spectra
spc.plot_spectrum(0,0,plot_fit=True)

# Show an image of the best-fit velocity
spc.mapplot.plane = spc.parcube[2,:,:]
spc.mapplot(estimator=None)

# running in script mode, the figures won't show by default on some systems
import pylab as pl
pl.show()
