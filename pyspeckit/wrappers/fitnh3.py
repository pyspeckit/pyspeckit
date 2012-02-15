"""
==================
NH3 fitter wrapper
==================

Wrapper to fit ammonia spectra.  Generates a reasonable guess at the position and velocity using a gaussian fit
"""
import pyspeckit
from matplotlib import pyplot 
import copy
import random

title_dict = {'oneone':'NH$_3(1,1)$', 
'twotwo':'NH$_3(2,2)$', 
'threethree':'NH$_3(3,3)$', 
'fourfour':'NH$_3(4,4)$', }

def fitnh3tkin(input_dict, dobaseline=True, baselinekwargs={}, crop=False, guessline='twotwo',
        tex=15,tkin=20,column=15.0,fortho=0.66, tau=None, thin=False, quiet=False, doplot=True, fignum=1,
        guessfignum=2, smooth=False, scale_keyword=None, rebase=False, npeaks=1, guesses=None,
        **kwargs): 
    """
    Given a dictionary of filenames and lines, fit them together
    e.g. {'oneone':'G000.000+00.000_nh3_11.fits'}
    """
    spdict = dict([ (linename,pyspeckit.Spectrum(value, scale_keyword=scale_keyword)) 
        if type(value) is str else (linename,value) for linename, value in input_dict.iteritems() ])
    splist = spdict.values()

    for sp in splist: # required for plotting, cropping
        sp.xarr.convert_to_unit('km/s')

    if crop and len(crop) == 2:
        for sp in splist:
            sp.crop(*crop)

    if dobaseline:
        for sp in splist:
            sp.baseline(**baselinekwargs)

    if smooth and type(smooth) is int:
        for sp in splist:
            sp.smooth(smooth)

    spdict[guessline].specfit(fittype='gaussian', negamp=False, vheight=False)
    ampguess,vguess,widthguess = spdict[guessline].specfit.modelpars
    if widthguess < 0:
        raise ValueError("Width guess was < 0.  This is impossible.")
    print "RMS guess (errspec): ",spdict[guessline].specfit.errspec.mean()
    print "RMS guess (residuals): ",spdict[guessline].specfit.residuals.std()
    errguess = spdict[guessline].specfit.residuals.std()

    if rebase:
        # redo baseline subtraction excluding the centroid +/- about 20 km/s 
        vlow = spdict[guessline].specfit.modelpars[1]-(19.8+spdict[guessline].specfit.modelpars[2]*2.35)
        vhigh = spdict[guessline].specfit.modelpars[1]+(19.8+spdict[guessline].specfit.modelpars[2]*2.35)
        for sp in splist:
            sp.baseline(exclude=[vlow,vhigh], **baselinekwargs)
    
    for sp in splist:
        sp.error[:] = errguess

    spdict[guessline].plotter(figure=guessfignum)
    spdict[guessline].specfit.plot_fit()

    spectra = pyspeckit.Spectra(splist)
    spectra.specfit.npeaks = npeaks

    if tau is not None:
        if guesses is None:
            guesses = [a for i in xrange(npeaks) for a in (tkin+random.random()*i,tex,tau+random.random()*i,widthguess+random.random()*i,vguess+random.random()*i,fortho)]
        spectra.specfit(fittype='ammonia_tau',quiet=quiet,multifit=True,guesses=guesses, thin=thin, **kwargs)
    else:
        if guesses is None:
            guesses = [a for i in xrange(npeaks) for a in (tkin+random.random()*i,tex,column+random.random()*i,widthguess+random.random()*i,vguess+random.random()*i,fortho)]
        spectra.specfit(fittype='ammonia',quiet=quiet,multifit=True,guesses=guesses, thin=thin, **kwargs)

    if doplot:
        plot_nh3(spdict,spectra,fignum=fignum)

    return spdict,spectra

def plot_nh3(spdict,spectra,fignum=1, show_components=False, residfignum=None, **plotkwargs):
    """
    Plot the results from a multi-nh3 fit
    """ 
    spectra.plotter.figure = pyplot.figure(fignum)
    spectra.plotter.axis = spectra.plotter.figure.gca()
    pyplot.clf()
    splist = spdict.values()

    for sp in splist:
        sp.xarr.convert_to_unit('km/s',quiet=True)
        sp.specfit.fitter = copy.copy(spectra.specfit.fitter)
        sp.specfit.modelpars = spectra.specfit.modelpars
        sp.specfit.npeaks = spectra.specfit.npeaks
        sp.specfit.fitter.npeaks = spectra.specfit.npeaks
        if spectra.specfit.modelpars is not None:
            sp.specfit.model = sp.specfit.fitter.n_ammonia(pars=spectra.specfit.modelpars, parnames=spectra.specfit.fitter.parnames)(sp.xarr)

    if len(splist) == 2:
        axdict = { 'oneone':pyplot.subplot(211), 'twotwo':pyplot.subplot(212) }
    elif len(splist) == 3:
        axdict = { 'oneone':pyplot.subplot(211), 'twotwo':pyplot.subplot(223), 'threethree':pyplot.subplot(224), 'fourfour':pyplot.subplot(224) }
    elif len(splist) == 4:
        axdict = { 'oneone':pyplot.subplot(221), 'twotwo':pyplot.subplot(222), 'threethree':pyplot.subplot(223), 'fourfour':pyplot.subplot(224) }
    for linename,sp in spdict.iteritems():
        sp.plotter.axis=axdict[linename] # permanent
        sp.plotter(axis=axdict[linename],title=title_dict[linename], **plotkwargs)
        sp.specfit.Spectrum.plotter = sp.plotter
        sp.specfit.selectregion(reset=True)
        if sp.specfit.modelpars is not None:
            sp.specfit.plot_fit(annotate=False, show_components=show_components)
    if spdict['oneone'].specfit.modelpars is not None:
        spdict['oneone'].specfit.annotate(labelspacing=0.05,prop={'size':'small','stretch':'extra-condensed'},frameon=False)

    if residfignum is not None:
        pyplot.figure(residfignum)
        pyplot.clf()
        if len(splist) == 2:
            axdict = { 'oneone':pyplot.subplot(211), 'twotwo':pyplot.subplot(212) }
        elif len(splist) == 3:
            axdict = { 'oneone':pyplot.subplot(211), 'twotwo':pyplot.subplot(223), 'threethree':pyplot.subplot(224), 'fourfour':pyplot.subplot(224) }
        elif len(splist) == 4:
            axdict = { 'oneone':pyplot.subplot(221), 'twotwo':pyplot.subplot(222), 'threethree':pyplot.subplot(223), 'fourfour':pyplot.subplot(224) }
        for linename,sp in spdict.iteritems():
            sp.specfit.plotresiduals(axis=axdict[linename])



def fitnh3(spectrum, vrange=[-100,100], vrangeunits='km/s', quiet=False,
        Tex=20,Tkin=15,column=1e15,fortho=1.0, tau=None): 

    if vrange:
        spectrum.xarr.convert_to_unit(vrangeunits)
        spectrum.crop(*vrange, units=vrangeunits)

    spectrum.specfit(fittype='gaussian',negamp=False)
    ampguess,vguess,widthguess = spectrum.specfit.modelpars

    if tau is None:
        spectrum.specfit(fittype='ammonia',quiet=quiet,multifit=True,guesses=[Tex,Tkin,column,widthguess,vguess,fortho])
    else:
        spectrum.specfit(fittype='ammonia_tau',quiet=quiet,multifit=True,guesses=[Tex,Tkin,tau,widthguess,vguess,fortho])

    return spectrum


def BigSpectrum_to_NH3dict(sp, vrange=None):
    """
    A rather complicated way to make the spdicts above given a spectrum...
    """

    spdict = {}
    for linename,freq in pyspeckit.spectrum.models.ammonia.freq_dict.iteritems():
        if vrange is not None:
            freq_test_low  = freq - freq * vrange[0]/pyspeckit.units.speedoflight_kms
            freq_test_high = freq - freq * vrange[1]/pyspeckit.units.speedoflight_kms
        else:
            freq_test_low = freq_test_high = freq

        if (sp.xarr.as_unit('Hz').in_range(freq_test_low) or
                sp.xarr.as_unit('Hz').in_range(freq_test_high)):
            spdict[linename] = sp.copy()
            spdict[linename].xarr.convert_to_unit('GHz')
            spdict[linename].xarr.refX = freq
            spdict[linename].xarr.refX_units = 'Hz'
            #spdict[linename].baseline = copy.copy(sp.baseline)
            #spdict[linename].baseline.Spectrum = spdict[linename]
            #spdict[linename].specfit = copy.copy(sp.specfit)
            #spdict[linename].specfit.fitter = copy.copy(sp.specfit.fitter)
            #spdict[linename].specfit.Spectrum = spdict[linename]
            #spdict[linename].specfit._full_model()
            spdict[linename].xarr.convert_to_unit('km/s')
            if vrange is not None:
                try:
                    spdict[linename].crop(*vrange, units='km/s')
                except IndexError:
                    # if the freq in range, but there's no data in range, remove
                    spdict.pop(linename)

    return spdict
                    
def plotter_override(sp, vrange=None, **kwargs):
    """
    Do plot_nh3 with syntax similar to plotter()
    """

    spdict = BigSpectrum_to_NH3dict(sp, vrange=vrange)

    if len(spdict) not in (2,3,4):
        raise ValueError("Not enough lines; don't need to use the NH3 plot wrapper")

    plot_nh3(spdict, sp, **kwargs)
