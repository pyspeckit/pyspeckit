import numpy as np
import pyspeckit
from astropy import units as u
from pyspeckit.spectrum.models import ammonia_constants, ammonia, ammonia_hf
from pyspeckit.spectrum.models.ammonia_constants import freq_dict
from pyspeckit.spectrum.units import SpectroscopicAxis, SpectroscopicAxes

# Step 1. Generate a synthetic spectrum.  Already have a real spectrum?  Skip
# to step 2!
# Generate a synthetic spectrum based off of 3 NH3 lines
# Note that they are converted to GHz first
xarr11 = SpectroscopicAxis(np.linspace(-30, 30, 100)*u.km/u.s,
                           velocity_convention='radio',
                           refX=freq_dict['oneone']).as_unit(u.GHz)
xarr22 = SpectroscopicAxis(np.linspace(-40, 40, 100)*u.km/u.s,
                           velocity_convention='radio',
                           refX=freq_dict['twotwo']).as_unit(u.GHz)
xarr33 = SpectroscopicAxis(np.linspace(-50, 50, 100)*u.km/u.s,
                           velocity_convention='radio',
                           refX=freq_dict['threethree']).as_unit(u.GHz)
# Merge the three X-axes into a single axis
xarr = SpectroscopicAxes([xarr11,xarr22,xarr33])

# Compute a synthetic model that is made of two temperature components with
# identical velocities
synthspec = (ammonia.ammonia(xarr, trot=20, ntot=15, fortho=0.5, xoff_v=0.0,
                             width=1.0) +
             ammonia.ammonia(xarr, trot=50, ntot=14, fortho=0.5, xoff_v=0.0,
                             width=1.0))

# Create the Spectrum object
spectrum = pyspeckit.Spectrum(xarr=xarr, data=synthspec, header={})


# Step 2.  You have a spectrum.
# plot it
spectrum.plotter()

# Use the multi-tex/multi-tau model generator to build up a model function
# You can use any set of oneone, twotwo, ..., eighteight (no 9-9 or higher)
# This sets the number of parameters to be fit: 2+2*(n_transitions)
fitter = ammonia_hf.nh3_vtau_multimodel_generator(['oneone', 'twotwo',
                                                   'threethree'])

# Register the fitter - i.e., tell pyspeckit where it is and how to use it
spectrum.specfit.Registry.add_fitter('nh3_vtau_123', fitter, fitter.npars)
# These are the parameter names, approximately:
# parnames=['center','width','Tex11','tau11','Tex22','tau22','Tex33','tau33'],

# Need to give some input guesses.  We start with something wrong-ish: -5 km/s,
# 1.2 km/s width, and 15 K + 0.5 tau for all 3 lines
guesses = [-5, 1.2, 15, 0.5, 15, 0.5, 15, 0.5,]

# Plot up the guessed model
spectrum.plotter.axis.plot(spectrum.xarr,
                           fitter.n_modelfunc(guesses)(spectrum.xarr), 'b')

# Run the fit!
spectrum.specfit(fittype='nh3_vtau_123', guesses=guesses)

# display the correct and fitted answers
print("Low column version:")
def printthings(ammonia=ammonia.ammonia, xarr=xarr):
    print("Real optical depths of component 1: ",[ammonia(xarr, trot=20,
                                                          ntot=15, fortho=0.5,
                                                          xoff_v=0.0,
                                                          width=1.0,
                                                          return_tau=True)[x]
                                                  for x in ['oneone', 'twotwo',
                                                            'threethree']])
    print("Real optical depths of component 2: ",[ammonia(xarr, trot=50,
                                                          ntot=14, fortho=0.5,
                                                          xoff_v=0.0,
                                                          width=1.0,
                                                          return_tau=True)[x]
                                                  for x in ['oneone', 'twotwo',
                                                            'threethree']])
printthings()

print("Fitted parameters: ",spectrum.specfit.parinfo)

# It works, but the covariances between tex and tau are large.
# So, another example with higher tau (and therefore... less degenerate?)

synthspec = (ammonia.ammonia(xarr, trot=20, ntot=16, fortho=0.5, xoff_v=0.0,
                             width=1.0) +
             ammonia.ammonia(xarr, trot=50, ntot=15, fortho=0.5, xoff_v=0.0,
                             width=1.0))
spectrum2 = pyspeckit.Spectrum(xarr=xarr, data=synthspec, header={})
spectrum2.plotter()
spectrum2.specfit.Registry.add_fitter('nh3_vtau_123', fitter, fitter.npars)
spectrum2.specfit(fittype='nh3_vtau_123', guesses=guesses)

# We can also examine what tau really should have been... kinda.
print("High column version:")
def printthings(ammonia=ammonia.ammonia, xarr=xarr):
    print("Real optical depths of component 1: ",[ammonia(xarr, trot=20,
                                                          ntot=16, fortho=0.5,
                                                          xoff_v=0.0,
                                                          width=1.0,
                                                          return_tau=True)[x]
                                                  for x in ['oneone', 'twotwo',
                                                            'threethree']])
    print("Real optical depths of component 2: ",[ammonia(xarr, trot=50,
                                                          ntot=15, fortho=0.5,
                                                          xoff_v=0.0,
                                                          width=1.0,
                                                          return_tau=True)[x]
                                                  for x in ['oneone', 'twotwo',
                                                            'threethree']])
printthings()
print("Fitted parameters: ",spectrum2.specfit.parinfo)
