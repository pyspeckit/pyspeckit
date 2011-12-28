"""
====================
N2H+ fitter wrapper
====================

Wrapper to fit N2H+ using RADEX models.  This is meant to be used from the command line, e.g.:

    python n2hp_wrapper.py file.fits

and therefore has no independently defined functions.
"""
import pyspeckit
import numpy as np
import pyfits
from pyspeckit.spectrum import models
from pyspeckit.spectrum.models import n2hp

# EDIT THIS LINE:
path_to_radex = '/Users/adam/work/n2hp/'

# create the n2hp Radex fitter
# This step cannot be easily generalized: the user needs to read in their own grids
texgrid1 = pyfits.getdata(path_to_radex+'1-2_T=5to55_lvg_tex1.fits')
taugrid1 = pyfits.getdata(path_to_radex+'1-2_T=5to55_lvg_tau1.fits')
texgrid2 = pyfits.getdata(path_to_radex+'1-2_T=5to55_lvg_tex2.fits')
taugrid2 = pyfits.getdata(path_to_radex+'1-2_T=5to55_lvg_tau2.fits')
hdr    = pyfits.getheader(path_to_radex+'1-2_T=5to55_lvg_tau2.fits')

# this deserves a lot of explanation:
# models.n2hp.n2hp_radex is the MODEL that we are going to fit
# models.model.SpectralModel is a wrapper to deal with parinfo, multiple peaks,
# and annotations
# all of the parameters after the first are passed to the model function 
n2hp_radex_fitter = models.model.SpectralModel(
        models.n2hp.n2hp_radex, 4,
        parnames=['density','column','center','width'], 
        parvalues=[4,12,0,1],
        parlimited=[(True,True), (True,True), (False,False), (True,False)], 
        parlimits=[(1,8), (11,16), (0,0), (0,0)],
        parsteps=[0.01,0.01,0,0],
        fitunits='Hz',
        texgrid=((4,5,texgrid1),(14,15,texgrid2)), # specify the frequency range over which the grid is valid (in GHz)
        taugrid=((4,5,taugrid1),(14,15,taugrid2)),
        hdr=hdr,
        shortvarnames=("n","N","v","\\sigma"), # specify the parameter names (TeX is OK)
        grid_vwidth_scale=False,
        )


if __name__ == "__main__":

    import optparse
    parser=optparse.OptionParser()
    parser.add_option("--scalekeyword",help="Scale the data before fitting?",default='ETAMB')
    parser.add_option("--vmin",help="Mininum velocity to include",default=-50)
    parser.add_option("--vmax",help="Maximum velocity to include",default=50)
    parser.add_option("--vguess",help="Velocity guess",default=3.75)
    parser.add_option("--smooth",help="Smooth the spectrum",default=None)

    options,args = parser.parse_args()

    if len(args) == 1:
        fn = args[0]
    else:
        print "Couldn't load a file!  Crash time!"


    sp = pyspeckit.Spectrum(fn,scale_keyword=options.scalekeyword)

    sp.xarr.convert_to_unit('km/s')
    sp.crop(options.vmin,options.vmax)
    if options.smooth: 
        sp.smooth(options.smooth) 

    sp.specfit() # determine errors
    sp.error = np.ones(sp.data.shape)*sp.specfit.residuals.std()
    sp.baseline(excludefit=True)

    sp.Registry.add_fitter('n2hp_radex',
            n2hp_radex_fitter,4,multisingle='multi')
    sp.Registry.add_fitter('n2hp_vtau',
            n2hp.n2hp_vtau_fitter,4,multisingle='multi')

    sp.plotter(figure=1)
    sp.specfit(fittype='n2hp_radex',multifit=True,guesses=[4,12,options.vguess,0.43],quiet=False)

    sp.plotter(figure=2)
    sp.specfit(fittype='n2hp_vtau',multifit=True,guesses=[15,2,options.vguess,0.43],quiet=False)
