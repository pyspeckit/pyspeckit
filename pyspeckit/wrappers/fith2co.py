"""
==================
H2CO fitter wrapper
==================

Wrapper to fit formaldehyde spectra.  
"""
import pyspeckit
from matplotlib import pyplot 
import copy
import random

title_dict = {'oneone':'H$_2$CO 1$_{11}$-1$_{10}$', 
'twotwo':'H$_2$CO 2$_{12}$-2$_{11}$',
'threethree':'H$_2$CO 3$_{23}$-3$_{22}$'}

def plot_h2co(spdict,spectra, fignum=1, show_components=False, residfignum=None, **plotkwargs):
    """
    Plot the results from a multi-h2co fit
    """ 
    spectra.plotter.figure = pyplot.figure(fignum)
    spectra.plotter.axis = spectra.plotter.figure.gca()
    pyplot.clf()
    splist = spdict.values()

    for sp in splist:
        sp.xarr.convert_to_unit('km/s',quiet=True)
        if hasattr(spectra.specfit,'fitter'):
            sp.specfit.fitter = copy.copy(spectra.specfit.fitter)
            sp.specfit.modelpars = spectra.specfit.modelpars
            sp.specfit.npeaks = spectra.specfit.npeaks
            sp.specfit.fitter.npeaks = spectra.specfit.npeaks
            if spectra.specfit.modelpars is not None:
                sp.specfit.model = sp.specfit.fitter.n_modelfunc(pars=spectra.specfit.modelpars, **spectra.specfit.fitter.modelfunc_kwargs)(sp.xarr)

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


def BigSpectrum_to_H2COdict(sp, vrange=None):
    """
    A rather complicated way to make the spdicts above given a spectrum...
    """

    spdict = {}
    for linename,freq in pyspeckit.spectrum.models.formaldehyde.central_freq_dict.iteritems():
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
            spdict[linename].specfit = sp.specfit.copy(parent=spdict[linename])
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
    Do plot_h2co with syntax similar to plotter()
    """

    spdict = BigSpectrum_to_H2COdict(sp, vrange=vrange)

    if len(spdict) not in (2,3,4):
        raise ValueError("Not enough lines; don't need to use the H2CO plot wrapper")

    plot_h2co(spdict, sp, **kwargs)
