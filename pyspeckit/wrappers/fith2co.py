"""
===================
H2CO fitter wrapper
===================

Wrapper to fit formaldehyde spectra.
"""
from __future__ import print_function
from .. import spectrum
from ..spectrum import units
import copy
from astropy import units as u
from six import iteritems

title_dict = {'oneone':'H$_2$CO 1$_{11}$-1$_{10}$',
              'twotwo':'H$_2$CO 2$_{12}$-2$_{11}$',
              'threethree':'H$_2$CO 3$_{23}$-3$_{22}$'
              }

def plot_h2co(spdict, spectra, fignum=1, show_components=False,
              show_hyperfine_components=False, residfignum=None, annotate=None,
              clear=True, residkwargs={}, plot_fit_kwargs={}, residclear=True,
              resid_overlay=False, resid_yoffsets=None,
              **plotkwargs):
    """
    Plot the results from a multi-h2co fit
    """
    from matplotlib import pyplot
    spectra.plotter.figure = pyplot.figure(fignum)
    spectra.plotter.axis = spectra.plotter.figure.gca()
    if clear:
        spectra.plotter.figure.clf()
    splist = spdict.values()

    for sp in splist:
        sp.xarr.convert_to_unit('km/s',quiet=True)
        if hasattr(spectra.specfit,'fitter'):
            sp.specfit.fitter = copy.copy(spectra.specfit.fitter)
            sp.specfit.modelpars = spectra.specfit.modelpars
            sp.specfit.npeaks = spectra.specfit.npeaks
            sp.specfit.fitter.npeaks = spectra.specfit.npeaks
            if spectra.specfit.modelpars is not None:
                mf = sp.specfit.fitter.n_modelfunc
                kw = spectra.specfit.fitter.modelfunc_kwargs
                sp.specfit.model = mf(pars=spectra.specfit.modelpars,
                                      **kw)(sp.xarr)

    if len(splist) == 2:
        axdict = {'oneone':pyplot.subplot(211),
                  'twotwo':pyplot.subplot(212)}
    elif len(splist) == 3:
        axdict = {'oneone':pyplot.subplot(211),
                  'twotwo':pyplot.subplot(223),
                  'threethree':pyplot.subplot(224)}
    elif len(splist) == 4:
        axdict = {'oneone':pyplot.subplot(221),
                  'twotwo':pyplot.subplot(222),
                  'threethree':pyplot.subplot(223),
                  'fourfour':pyplot.subplot(224)}
    for linename,sp in iteritems(spdict):
        sp.plotter.axis=axdict[linename] # permanent
        sp.plotter(axis=axdict[linename],
                   title=title_dict[linename],
                   clear=clear,
                   **plotkwargs)
        sp.specfit.Spectrum.plotter = sp.plotter
        #sp.specfit.selectregion(reset=True)
        if sp.specfit.modelpars is not None:
            sp.specfit.plot_fit(annotate=False,
                                show_components=show_components,
                                show_hyperfine_components=show_hyperfine_components,
                                **plot_fit_kwargs)
        sp.plotter.reset_limits()
    if spdict['oneone'].specfit.modelpars is not None and annotate:
        spdict['oneone'].specfit.annotate(labelspacing=0.05,
                                          prop={'size':'small',
                                                'stretch':'extra-condensed'},
                                          frameon=False)

    residaxdict = None
    if residfignum is not None:
        pyplot.figure(residfignum)
        if residclear:
            pyplot.clf()
        if len(splist) == 2:
            residaxdict = {'oneone':pyplot.subplot(211),
                           'twotwo':pyplot.subplot(212)}
        elif len(splist) == 3:
            residaxdict = {'oneone':pyplot.subplot(211),
                           'twotwo':pyplot.subplot(223),
                           'threethree':pyplot.subplot(224),
                           'fourfour':pyplot.subplot(224)}
        elif len(splist) == 4:
            residaxdict = {'oneone':pyplot.subplot(221),
                           'twotwo':pyplot.subplot(222),
                           'threethree':pyplot.subplot(223),
                           'fourfour':pyplot.subplot(224)}
    elif resid_overlay:
        residaxdict = axdict
        residclear = False # override defaults...
        residfignum = fignum

    if residaxdict is not None:
        for linename,sp in iteritems(spdict):
            sp.specfit.Spectrum.plotter = sp.plotter
            try:
                yoffset = resid_yoffsets[linename]
            except TypeError:
                yoffset = 0.0
            sp.specfit.plotresiduals(axis=residaxdict[linename],
                                     figure=residfignum,
                                     clear=residclear,
                                     set_limits=False,
                                     label=False,
                                     yoffset=yoffset,
                                     **residkwargs)
        spectra.residaxdict = residaxdict

    spectra.axisdict = axdict
    spectra.plotter.axis = axdict['oneone']
    spectra.specfit.fitleg = spdict['oneone'].specfit.fitleg


def BigSpectrum_to_H2COdict(sp, vrange=None):
    """
    A rather complicated way to make the spdicts above given a spectrum...
    """

    spdict = {}
    for linename,freq in iteritems(spectrum.models.formaldehyde.central_freq_dict):
        if vrange is not None:
            freq_test_low  = freq - freq * vrange[0]/units.speedoflight_kms
            freq_test_high = freq - freq * vrange[1]/units.speedoflight_kms
        else:
            freq_test_low = freq_test_high = freq

        if (sp.xarr.as_unit('Hz').in_range(freq_test_low*u.Hz) or
                sp.xarr.as_unit('Hz').in_range(freq_test_high*u.Hz)):
            spdict[linename] = sp.copy(deep=True)
            spdict[linename].xarr.convert_to_unit('GHz')
            spdict[linename].xarr.refX = freq
            spdict[linename].xarr.refX_unit = 'Hz'
            #spdict[linename].baseline = copy.copy(sp.baseline)
            #spdict[linename].baseline.Spectrum = spdict[linename]
            spdict[linename].specfit = sp.specfit.copy(parent=spdict[linename])
            spdict[linename].xarr.convert_to_unit('km/s')
            if vrange is not None:
                try:
                    spdict[linename].crop(*vrange, unit='km/s')
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

    return spdict
