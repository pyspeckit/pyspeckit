import pyspeckit
import numpy as np
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
from pyspeckit.spectrum import models

# create the Formaldehyde Radex fitter
# This step cannot be easily generalized: the user needs to read in their own grids
texgrid1 = pyfits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_321_220_temperature_para_tex1.fits')
taugrid1 = pyfits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_321_220_temperature_para_tau1.fits')
texgrid2 = pyfits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_321_220_temperature_para_tex2.fits')
taugrid2 = pyfits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_321_220_temperature_para_tau2.fits')
hdr = pyfits.getheader('/Users/adam/work/h2co/radex/thermom/303-202_321_220_temperature_para_tau2.fits')

texgrid1b = pyfits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_322-221_temperature_para_tex1.fits')
taugrid1b = pyfits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_322-221_temperature_para_tau1.fits')
texgrid2b = pyfits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_322-221_temperature_para_tex2.fits')
taugrid2b = pyfits.getdata('/Users/adam/work/h2co/radex/thermom/303-202_322-221_temperature_para_tau2.fits')
hdrb = pyfits.getheader('/Users/adam/work/h2co/radex/thermom/303-202_322-221_temperature_para_tau2.fits')

# this deserves a lot of explanation:
# models.formaldehyde.formaldehyde_radex is the MODEL that we are going to fit
# models.model.SpectralModel is a wrapper to deal with parinfo, multiple peaks,
# and annotations
# all of the parameters after the first are passed to the model function 
formaldehyde_radex_fitter_b = models.model.SpectralModel(
        models.formaldehyde_mm.formaldehyde_mm_radex, 5,
        parnames=['temperature','column','density','center','width'], 
        parvalues=[50,12,4.5,0,1],
        parlimited=[(True,True), (True,True), (True,True), (False,False), (True,False)], 
        parlimits=[(5,205), (10,17), (2,7), (0,0), (0,0)],
        parsteps=[0.01,0.01,0.1,0,0],
        fitunits='Hz',
        texgrid=((218.2,218.3,texgrid1b),(218.4,218.55,texgrid2b)), # specify the frequency range over which the grid is valid (in GHz)
        taugrid=((218.2,218.3,taugrid1b),(218.4,218.55,taugrid2b)),
        hdr=hdrb,
        shortvarnames=("T","N","n","v","\\sigma"), # specify the parameter names (TeX is OK)
        grid_vwidth=5.0,
        )

formaldehyde_radex_fitter = models.model.SpectralModel(
        models.formaldehyde_mm.formaldehyde_mm_radex, 5,
        parnames=['temperature','column','density','center','width'], 
        parvalues=[50,12,4.5,0,1],
        parlimited=[(True,True), (True,True), (True,True), (False,False), (True,False)], 
        parlimits=[(5,205), (10,17), (2,7), (0,0), (0,0)],
        parsteps=[0.01,0.01,0.1,0,0],
        fitunits='Hz',
        texgrid=((218.2,218.3,texgrid1),(218.7,218.8,texgrid2)), # specify the frequency range over which the grid is valid (in GHz)
        taugrid=((218.2,218.3,taugrid1),(218.7,218.8,taugrid2)),
        hdr=hdr,
        shortvarnames=("T","N","n","v","\\sigma"), # specify the parameter names (TeX is OK)
        grid_vwidth=5.0,
        )

formaldehyde_radex_fitter_both = models.model.SpectralModel(
        models.formaldehyde_mm.formaldehyde_mm_radex, 5,
        parnames=['temperature','column','density','center','width'], 
        parvalues=[50,12,4.5,0,1],
        parlimited=[(True,True), (True,True), (True,True), (False,False), (True,False)], 
        parlimits=[(5,205), (10,17), (2,7), (0,0), (0,0)],
        parsteps=[0.01,0.01,0.1,0,0],
        fitunits='Hz',
        texgrid=((218.2,218.3,texgrid1b),(218.4,218.55,texgrid2b),(218.7,218.8,texgrid2)), # specify the frequency range over which the grid is valid (in GHz)
        taugrid=((218.2,218.3,taugrid1b),(218.4,218.55,taugrid2b),(218.7,218.8,taugrid2)),
        hdr=hdrb,
        shortvarnames=("T","N","n","v","\\sigma"), # specify the parameter names (TeX is OK)
        grid_vwidth=5.0,
        )


if __name__ == "__main__":

    import optparse
    parser=optparse.OptionParser()
    parser.add_option("--scalekeyword",help="Scale the data before fitting?",default='ETAMB')
    parser.add_option("--vguess",help="Velocity guess",default=0.0)
    parser.add_option("--smooth6cm",help="Smooth the 6cm spectrum",default=3)
    parser.add_option("--pymc",help="pymc?",default=False)

    options,args = parser.parse_args()

    if len(args) == 1:
        sp = pyspeckit.Spectrum(args[0])
    else:
        sp = pyspeckit.readers.read_class.class_to_spectra('spectra_off_30_-60.apex',apex=True)
        sp.error[:] = sp.stats((2.183e2,2.184e2))['std']

    sp.Registry.add_fitter('formaldehyde_mm_radex',
            formaldehyde_radex_fitter,5,multisingle='multi',            )
    sp.Registry.add_fitter('formaldehyde_mm_radex_b',
            formaldehyde_radex_fitter_b,5,multisingle='multi',            )
    sp.Registry.add_fitter('formaldehyde_mm_radex_both',
            formaldehyde_radex_fitter_both,5,multisingle='multi',            )

    sp1 = sp.copy()
    sp2 = sp.copy()
    sp3 = sp.copy()
    sp1.plotter(figure=1)
    sp1.specfit(fittype='formaldehyde_mm_radex',multifit=True,guesses=[100,13.2,4.5,options.vguess,7.0],
            limits=[(20,200),(11,15),(3,5.5),(-5,5),(2,15)], limited=[(True,True)]*5, 
            fixed=[False,False,True,False,False], quiet=False,)
    sp2.plotter(figure=2)
    sp2.specfit(fittype='formaldehyde_mm_radex_b',multifit=True,guesses=[100,13.2,4.5,options.vguess,7.0],
            limits=[(20,200),(11,15),(3,5.5),(-5,5),(2,15)], limited=[(True,True)]*5, 
            fixed=[False,False,True,False,False], quiet=False,)
    sp3.plotter(figure=3)
    sp3.specfit(fittype='formaldehyde_mm_radex_both',multifit=True,guesses=[100,13.2,4.5,options.vguess,7.0],
            limits=[(20,200),(11,15),(3,5.5),(-5,5),(2,15)], limited=[(True,True)]*5, 
            fixed=[False,False,False,False,False], quiet=False,)
    sp3.specfit.parinfo.TEMPERATURE0.error=50
    sp3.specfit.parinfo.COLUMN0.error=2
    sp3.specfit.parinfo.DENSITY0.error=2
    sp3.specfit.parinfo.WIDTH0.fixed=True
    sp3.specfit.parinfo.CENTER0.fixed=True

    if options.pymc:
        spmc = sp3.specfit.get_pymc(use_fitted_values=True)
        spmc.sample(10000)
        import agpy
        agpy.pymc_plotting.hist2d(spmc, 'TEMPERATURE0', 'DENSITY0', doerrellipse=False, clear=True, fignum=4)

