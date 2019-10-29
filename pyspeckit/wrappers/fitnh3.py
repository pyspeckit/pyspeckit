"""
NH3 fitter wrapper
==================

Wrapper to fit ammonia spectra.  Generates a reasonable guess at the position
and velocity using a gaussian fit

Example use:

.. code:: python

    import pyspeckit
    sp11 = pyspeckit.Spectrum('spec.nh3_11.dat', errorcol=999)
    sp22 = pyspeckit.Spectrum('spec.nh3_22.dat', errorcol=999)
    sp33 = pyspeckit.Spectrum('spec.nh3_33.dat', errorcol=999)
    sp11.xarr.refX = pyspeckit.spectrum.models.ammonia.freq_dict['oneone']
    sp22.xarr.refX = pyspeckit.spectrum.models.ammonia.freq_dict['twotwo']
    sp33.xarr.refX = pyspeckit.spectrum.models.ammonia.freq_dict['threethree']
    input_dict={'oneone':sp11, 'twotwo':sp22, 'threethree':sp33}
    spf = pyspeckit.wrappers.fitnh3.fitnh3tkin(input_dict)


Note that if you want to use the plotter wrapper with cubes, you need to do
something like the following, where the ``plot_special`` method of the stacked
``cubes`` object is set to the ``plotter_override`` function defined in the
fitnh3_wrapper code:

.. code:: python

    cubes.plot_special = pyspeckit.wrappers.fitnh3.plotter_override
    cubes.plot_special_kwargs = {'fignum':3, 'vrange':[55,135]}
    cubes.plot_spectrum(160,99)


"""
from __future__ import print_function
import warnings
from six.moves import xrange
from six import iteritems
import pyspeckit
from .. import spectrum
from ..spectrum.classes import Spectrum, Spectra
from ..spectrum import units
from ..spectrum.models import ammonia_constants
import numpy as np
import copy
import random
from astropy import log
from astropy import units as u

pyspeckit.spectrum.fitters.default_Registry.add_fitter('ammonia_tau_thin',
                                                       pyspeckit.spectrum.models.ammonia.ammonia_model_vtau_thin(),
                                                       5)

title_dict = {'oneone':'NH$_3(1, 1)$', 'twotwo':'NH$_3(2, 2)$',
              'threethree':'NH$_3(3, 3)$', 'fourfour':'NH$_3(4, 4)$',
              'fivefive':'NH$_3(5, 5)$', 'sixsix':'NH$_3(6, 6)$',
              'sevenseven':'NH$_3(7, 7)$', 'eighteight':'NH$_3(8, 8)$',
             }

def fitnh3tkin(input_dict, dobaseline=True, baselinekwargs={}, crop=False,
               cropunit=None, guessline='twotwo', tex=15, trot=20, column=15.0,
               fortho=0.66, tau=None, thin=False, quiet=False, doplot=True,
               fignum=1, guessfignum=2, smooth=False, scale_keyword=None,
               rebase=False, tkin=None, npeaks=1, guesses=None,
               fittype='ammonia',
               guess_error=True, plotter_wrapper_kwargs={}, **kwargs):
    """
    Given a dictionary of filenames and lines, fit them together
    e.g. {'oneone':'G000.000+00.000_nh3_11.fits'}

    Parameters
    ----------
    input_dict : dict
        A dictionary in which the keys are the ammonia line names (e.g.,
        'oneone', 'twotwo', etc) and the values are either Spectrum objects
        or filenames of spectra
    dobaseline : bool
        Fit and subtract a baseline prior to fitting the model?
        Keyword arguments to `pyspeckit.spectrum.Spectrum.baseline` are
        specified in ``baselinekwargs``.
    baselinekwargs : dict
        The keyword arguments for the baseline
    crop : bool or tuple
        A range of values to crop the spectrum to.  The units are specified by
        ``cropunit``; the default ``None`` will use pixels.  If False, no
        cropping will be performed.
    cropunit : None or astropy unit
        The unit for the crop parameter
    guess_error : bool
        Use the guess line to estimate the error in all spectra?
    plotter_wrapper_kwargs : dict
        Keyword arguments to pass to the plotter
    fittype: 'ammonia' or 'cold_ammonia'
        The fitter model to use.  This is overridden if `tau` is specified,
        in which case one of the `ammonia_tau` models is used (see source code)
    """
    if tkin is not None:
        if trot == 20 or trot is None:
            trot = tkin
        else:
            raise ValueError("Please specify trot, not tkin")
        warnings.warn("Keyword 'tkin' is deprecated; use trot instead", DeprecationWarning)

    spdict = dict([(linename, Spectrum(value, scale_keyword=scale_keyword))
                   if type(value) is str else (linename, value)
                   for linename, value in iteritems(input_dict)
                  ])
    splist = spdict.values()

    for transition, sp in spdict.items(): # required for plotting, cropping
        sp.xarr.convert_to_unit('km/s', velocity_convention='radio',
                                refX=pyspeckit.spectrum.models.ammonia.freq_dict[transition]*u.Hz,
                                quiet=True)

    if crop and len(crop) == 2:
        for sp in splist:
            sp.crop(*crop, unit=cropunit)

    if dobaseline:
        for sp in splist:
            sp.baseline(**baselinekwargs)

    if smooth and type(smooth) is int:
        for sp in splist:
            sp.smooth(smooth)

    spdict[guessline].specfit(fittype='gaussian', negamp=False, vheight=False,
                              guesses='moments')
    ampguess, vguess, widthguess = spdict[guessline].specfit.modelpars
    if widthguess < 0:
        raise ValueError("Width guess was < 0.  This is impossible.")
    print("RMS guess (errspec): ", spdict[guessline].specfit.errspec.mean())
    print("RMS guess (residuals): ", spdict[guessline].specfit.residuals.std())
    errguess = spdict[guessline].specfit.residuals.std()

    if rebase:
        # redo baseline subtraction excluding the centroid +/- about 20 km/s
        vlow = spdict[guessline].specfit.modelpars[1]-(19.8+spdict[guessline].specfit.modelpars[2]*2.35)
        vhigh = spdict[guessline].specfit.modelpars[1]+(19.8+spdict[guessline].specfit.modelpars[2]*2.35)
        for sp in splist:
            sp.baseline(exclude=[vlow, vhigh], **baselinekwargs)

    for sp in splist:
        if guess_error:
            sp.error[:] = errguess
        sp.xarr.convert_to_unit(u.GHz)

    if doplot:
        spdict[guessline].plotter(figure=guessfignum)
        spdict[guessline].specfit.plot_fit()

    spectra = Spectra(splist)
    spectra.specfit.npeaks = npeaks

    if tau is not None:
        if guesses is None:
            guesses = [a for i in xrange(npeaks) for a in
                       (trot+random.random()*i, tex, tau+random.random()*i,
                        widthguess+random.random()*i, vguess+random.random()*i,
                        fortho)]
        fittype = 'ammonia_tau_thin' if thin else 'ammonia_tau'
        spectra.specfit(fittype=fittype, quiet=quiet, guesses=guesses,
                        **kwargs)
    else:
        if guesses is None:
            guesses = [a for i in xrange(npeaks) for a in
                       (trot+random.random()*i, tex, column+random.random()*i,
                        widthguess+random.random()*i, vguess+random.random()*i,
                        fortho)]
        if thin:
            raise ValueError("'thin' keyword not supported for the generic ammonia model")
        spectra.specfit(fittype=fittype, quiet=quiet, guesses=guesses,
                        **kwargs)

    if doplot:
        plot_nh3(spdict, spectra, fignum=fignum, **plotter_wrapper_kwargs)

    return spdict, spectra

def plot_nh3(spdict, spectra, fignum=1, show_components=False,
             residfignum=None, show_hyperfine_components=True, annotate=True,
             axdict=None, figure=None,
             **plotkwargs):
    """
    Plot the results from a multi-nh3 fit

    spdict needs to be dictionary with form:
        'oneone': spectrum,
        'twotwo': spectrum,
        etc.
    """
    from matplotlib import pyplot
    if figure is None:
        spectra.plotter.figure = pyplot.figure(fignum)
        spectra.plotter.axis = spectra.plotter.figure.gca()

    splist = spdict.values()

    for transition, sp in spdict.items():
        sp.xarr.convert_to_unit('km/s', velocity_convention='radio',
                                refX=pyspeckit.spectrum.models.ammonia.freq_dict[transition]*u.Hz,
                                quiet=True)
        try:
            sp.specfit.fitter = copy.copy(spectra.specfit.fitter)
            sp.specfit.fitter.npeaks = spectra.specfit.npeaks
        except AttributeError:
            pass
        sp.specfit.modelpars = spectra.specfit.modelpars
        sp.specfit.parinfo = spectra.specfit.parinfo
        sp.specfit.npeaks = spectra.specfit.npeaks
        if spectra.specfit.modelpars is not None:
            sp.specfit.model = sp.specfit.fitter.n_ammonia(pars=spectra.specfit.modelpars, parnames=spectra.specfit.fitter.parnames)(sp.xarr)

    if axdict is None:
        axdict = make_axdict(splist, spdict)

    for linename, sp in iteritems(spdict):
        if linename not in axdict:
            raise NotImplementedError("Plot windows for {0} cannot "
                                      "be automatically arranged (yet)."
                                      .format(linename))
        sp.plotter.axis=axdict[linename] # permanent
        sp.plotter(axis=axdict[linename], title=title_dict[linename], **plotkwargs)
        sp.specfit.Spectrum.plotter = sp.plotter
        sp.specfit.selectregion(reset=True)
        if sp.specfit.modelpars is not None:
            sp.specfit.plot_fit(annotate=False, show_components=show_components,
                                show_hyperfine_components=show_hyperfine_components)
    if spdict['oneone'].specfit.modelpars is not None and annotate:
        spdict['oneone'].specfit.annotate(labelspacing=0.05,
                                          prop={'size':'small',
                                                'stretch':'extra-condensed'},
                                          frameon=False)

    if residfignum is not None:
        pyplot.figure(residfignum)
        pyplot.clf()
        axdict = make_axdict(splist, spdict)
        for linename, sp in iteritems(spdict):
            sp.specfit.plotresiduals(axis=axdict[linename])


def make_axdict(splist, spdict):
    from matplotlib import pyplot
    axdict = {}
    if len(splist) == 2:
        ii = 1
        for linename in ammonia_constants.line_names:
            if linename in spdict:
                axdict[linename] = pyplot.subplot(2,1,ii)
                ii+=1
    elif len(splist) == 3:
        ii = 1
        for linename in ammonia_constants.line_names:
            if linename in spdict:
                if ii == 1:
                    axdict[linename] = pyplot.subplot(2,1,ii)
                    ii+=2
                else:
                    axdict[linename] = pyplot.subplot(2,2,ii)
                    ii+=1
    elif len(splist) == 4:
        ii = 1
        for linename in ammonia_constants.line_names:
            if linename in spdict:
                axdict[linename] = pyplot.subplot(2,2,ii)
                ii+=1
    else:
        raise NotImplementedError("Plots with {0} subplots are not yet "
                                  "implemented.  Pull requests are "
                                  "welcome!".format(len(splist)))

    return axdict


def fitnh3(spectrum, vrange=[-100, 100], vrangeunit='km/s', quiet=False, Tex=20,
           trot=15, column=1e15, fortho=1.0, tau=None, Tkin=None,
           fittype='ammonia',
           spec_convert_kwargs={}):

    if Tkin is not None:
        if trot == 20 or trot is None:
            trot = Tkin
        else:
            raise ValueError("Please specify trot, not Tkin")
        warnings.warn("Keyword 'Tkin' is deprecated; use trot instead", DeprecationWarning)

    if vrange:
        spectrum.xarr.convert_to_unit(vrangeunit, **spec_convert_kwargs)
        spectrum.crop(*vrange, unit=vrangeunit)

    spectrum.specfit(fittype='gaussian', negamp=False, guesses='moments')
    ampguess, vguess, widthguess = spectrum.specfit.modelpars

    if tau is None:
        spectrum.specfit(fittype=fittype, quiet=quiet,
                         guesses=[Tex, trot, column, widthguess, vguess,
                                  fortho])
    else:
        spectrum.specfit(fittype='ammonia_tau', quiet=quiet,
                         guesses=[Tex, trot, tau, widthguess, vguess, fortho])

    return spectrum


def BigSpectrum_to_NH3dict(sp, vrange=None):
    """
    A rather complicated way to make the spdicts above given a spectrum...
    """

    sp.xarr.convert_to_unit('GHz')

    spdict = {}
    for linename, freq in iteritems(spectrum.models.ammonia.freq_dict):
        if not hasattr(freq, 'unit'):
            freq = freq*u.Hz
        if vrange is not None:
            freq_test_low  = freq - freq * vrange[0]/units.speedoflight_kms
            freq_test_high = freq - freq * vrange[1]/units.speedoflight_kms
        else:
            freq_test_low = freq_test_high = freq

        log.debug("line {2}: freq test low, high: {0}, {1}"
                  .format(freq_test_low, freq_test_high, linename))
        if (sp.xarr.as_unit('Hz').in_range(freq_test_low) or
                sp.xarr.as_unit('Hz').in_range(freq_test_high)):
            spdict[linename] = sp.copy(deep=True)
            spdict[linename].xarr.convert_to_unit('GHz')
            assert np.all(np.array(spdict[linename].xarr == sp.xarr,
                          dtype='bool'))
            spdict[linename].xarr.refX = freq
            spdict[linename].xarr.convert_to_unit('km/s',
                                                  velocity_convention='radio',
                                                  refX=pyspeckit.spectrum.models.ammonia.freq_dict[linename]*u.Hz,
                                                  quiet=True)
            np.testing.assert_array_almost_equal(spdict[linename].xarr.as_unit('GHz').value,
                                                 sp.xarr.value)
            log.debug("Line {0}={2}: {1}".format(linename, spdict[linename],
                                                 freq))
            if vrange is not None:
                try:
                    spdict[linename] = spdict[linename].slice(start=vrange[0],
                                                              stop=vrange[1],
                                                              unit='km/s')
                    log.debug("Successfully cropped {0} to {1}, freq = {2}, {3}"
                              .format(linename, vrange, freq,
                                      spdict[linename].xarr))
                    if len(spdict[linename]) == 0:
                        spdict.pop(linename)
                        log.debug("Removed {0} from spdict".format(linename))
                except IndexError:
                    # if the freq in range, but there's no data in range, remove
                    spdict.pop(linename)
        else:
            log.debug("Line {0} not in spectrum".format(linename))

    # this shouldn't be reachable, but there are reported cases where spdict
    # gets populated w/empty spectra, which leads to a failure in producing
    # their repr.  Since that on its own isn't a very helpful error message,
    # we'd rather return the bad spdict and see if the next function down the
    # line can survive with a questionable spdict...
    try:
        log.debug(str(spdict))
    except Exception as ex:
        log.debug(str(ex))

    return spdict

def plotter_override(sp, vrange=None, **kwargs):
    """
    Do plot_nh3 with syntax similar to plotter()
    """

    spdict = BigSpectrum_to_NH3dict(sp, vrange=vrange)
    log.debug("spdict: {0}".format(spdict))

    if len(spdict) > 4:
        raise ValueError("Too many lines ({0}) found.".format(len(spdict)))
    if len(spdict) not in (2, 3, 4):
        raise ValueError("Not enough lines; don't need to use the NH3 plot "
                         "wrapper.  If you think you are getting this message "
                         "incorrectly, check the velocity range (vrange "
                         "parameter) and make sure your spectrum overlaps with "
                         " it.")

    plot_nh3(spdict, sp, **kwargs)

    return spdict
