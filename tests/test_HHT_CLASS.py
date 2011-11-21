import pyspeckit
from pylab import *

if not 'interactive' in globals():
    interactive=False

sp = pyspeckit.Spectrum('10074-190_HCOp.fits')
sp.plotter(title=sp.header.get('OBJECT'))
sp.baseline(order=2)
# sp.xarr.refX=center_frequency=sp.header.get('RESTFREQ')
sp.xarr.convert_to_unit('km/s')
sp.plotter()
if interactive: raw_input('Wait - just plot it')
sp.plotter(xmin=-100,xmax=150)
sp.specfit(negamp=False,limitedmin=[True,True,False,True])
#print sp.specfit.guesses,sp.specfit.modelpars,sp.specfit.modelerrs,sp.specfit.errspec.mean(),sp.specfit.errspec.std()
if interactive: raw_input('Wait - show baseline before subtracting')
if sp.specfit.modelpars[0] > sp.specfit.modelerrs[0]:
    sp.baseline(excludefit=True,subtract=False,order=2,quiet=False)
    sp.specfit(negamp=False,limitedmin=[True,True,False,True])
    #print sp.specfit.guesses,sp.specfit.modelpars,sp.specfit.modelerrs,sp.specfit.errspec.mean(),sp.specfit.errspec.std()
    #print sp.specfit.fitter.annotations()
if interactive: raw_input('Wait - show baseline before subtracting')
sp.plotter()
sp.baseline(excludefit=True,subtract=True,order=2)
sp.specfit(negamp=False,limitedmin=[True,True,False,True])
if interactive: raw_input("Wait - we're done now")
#sp.baseline(excludefit=True,subtract=False,quiet=False,LoudDebug=True)
#raw_input('Wait')
