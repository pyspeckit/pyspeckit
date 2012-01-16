import pyspeckit
import matplotlib

if not 'interactive' in globals():
    interactive=False
if not 'savedir' in globals():
    savedir = ''

# load a FITS-compliant spectrum
spec = pyspeckit.Spectrum('10074-190_HCOp.fits')
# The units are originally frequency (check this by printing spec.xarr.units).
# I want to know the velocity.  Convert!
# Note that this only works because the reference frequency is set in the header
# this is no longer necessary!  #spec.xarr.frequency_to_velocity()
# Default conversion is to m/s, but we traditionally work in km/s
spec.xarr.convert_to_unit('km/s')
# plot it up!
spec.plotter()
# Subtract a baseline (the data is only 'mostly' reduced)
spec.baseline(interactive=True)
event1 = matplotlib.backend_bases.MouseEvent('button_press_event', spec.plotter.axis.figure.canvas,143,223,button=1)
event2 = matplotlib.backend_bases.MouseEvent('button_press_event', spec.plotter.axis.figure.canvas,400,223,button=1)
event3 = matplotlib.backend_bases.MouseEvent('button_press_event', spec.plotter.axis.figure.canvas,472,223,button=1)
event4 = matplotlib.backend_bases.MouseEvent('button_press_event', spec.plotter.axis.figure.canvas,787,223,button=1)
event5 = matplotlib.backend_bases.MouseEvent('button_press_event', spec.plotter.axis.figure.canvas,787,223,button=2)
spec.baseline.event_manager(event1)
spec.baseline.event_manager(event2)
spec.baseline.event_manager(event3)
spec.baseline.event_manager(event4)
spec.baseline.event_manager(event5)

spec.specfit(interactive=True)
event1 = matplotlib.backend_bases.KeyEvent('button_press_event', spec.plotter.axis.figure.canvas,x=463,y=570,key=2)
event2 = matplotlib.backend_bases.KeyEvent('button_press_event', spec.plotter.axis.figure.canvas,x=461,y=351,key=2)
event3 = matplotlib.backend_bases.KeyEvent('button_press_event', spec.plotter.axis.figure.canvas,x=403,y=256,key=1)
event4 = matplotlib.backend_bases.KeyEvent('button_press_event', spec.plotter.axis.figure.canvas,x=516,y=243,key=1)
event5 = matplotlib.backend_bases.KeyEvent('button_press_event', spec.plotter.axis.figure.canvas,x=597,y=257,key=3)
spec.specfit.event_manager(event1)
spec.specfit.event_manager(event2)
spec.specfit.event_manager(event3)
spec.specfit.event_manager(event4)
spec.specfit.event_manager(event5)

spec.baseline(excludefit=True)
spec.specfit(guesses=spec.specfit.modelpars)

spec.plotter.figure.savefig(savedir+"simple_fit_interactive_HCOp.png")
spec.specfit.residualaxis.figure.savefig(savedir+"simple_fit_interactive_HCOp_residuals.png")
