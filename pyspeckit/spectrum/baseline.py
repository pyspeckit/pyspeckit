from __future__ import print_function
import numpy as np
import matplotlib
from ..config import ConfigDescriptor as cfgdec
from .. import specwarnings
from .. import mpfit
from . import interactive
from . import history
from . import models
from astropy import log
import copy

interactive_help_message = """
(1) Left-click or press 1 (one) at two positions to select or add to the baseline fitting range - it will be
highlighted in green if the selection is successful.
    You can select regions to e/x/clude by pressing 'x' at two positions
(2) Middle or right click or press '2','m', '3', or 'd' to /d/isconnect and perform the fit.
    If you press '2','m', or middle-click, the baseline will be subtracted
    If you press '3','d', or right-click, the baseline will be plotted but not subtracted
"""

class Baseline(interactive.Interactive):
    """
    Class to measure and subtract baselines from spectra.

    While the term 'baseline' is generally used in the radio to refer to
    broad-band features in a spectrum not necessarily associated with a source,
    in this package it refers to general continuum fitting.  In principle,
    there's no reason to separate 'continuum' and 'spectral feature' fitting
    into different categories (both require some model, data, and optional
    weights when fitting).  In practice, however, 'continuum' is frequently
    something to be removed and ignored, while spectral features are the
    desired measurable quantity.  In order to accurately measure spectral
    features, it is necessary to allow baselines of varying complexity.

    The Baseline class has both interactive and command-based data selection
    features.  It can be used to fit both polynomial and power-law continua.
    Blackbody fitting is not yet implemented [12/21/2011].  Baseline fitting is
    a necessary prerequisite for Equivalent Width measurement.

    As you may observe in the comments on this code, this has been one of the
    buggiest and least adequately tested components of pyspeckit.  Bug reports
    are welcome.  (as of 1/15/2012, a major change has probably fixed most of
    the bugs, and the code base is much simpler)
    """
    def __init__(self, Spectrum):
        super(Baseline, self).__init__(Spectrum,
                interactive_help_message=interactive_help_message)
        self.baselinepars  = None
        self.order = None
        self.basespec = np.zeros(Spectrum.data.shape[0])
        #self.excludemask = np.zeros(Spectrum.data.shape[0],dtype='bool')
        self.OKmask = np.ones(Spectrum.data.shape[0],dtype='bool')
        #self.Spectrum = Spectrum
        self.Spectrum.plotter = Spectrum.plotter
        self.blleg = None
        self.click = 0
        #self.nclicks_b1 = 0
        #self.fitregion=[]
        #self.excludevelo = []
        #self.excludepix  = []
        self.subtracted = False
        #self.xmin = 0
        #self.xmax = Spectrum.data.shape[0] - 1
        #self.include = [self.xmin,self.xmax]
        #self.includevelo = [Spectrum.xarr[self.xmin],Spectrum.xarr[self.xmax]]
        self.powerlaw=False
        self.spline=False
        self._xfit_units = None
        self._plots = []

    @cfgdec
    def __call__(self, order=1, excludefit=False, save=True, exclude=None,
                 exclusionlevel=0.01, interactive=False, debug=False,
                 LoudDebug=False, fit_original=True,
                 baseline_fit_color='orange', clear_all_connections=True,
                 fit_plotted_area=True, highlight_fitregion=False,
                 reset_selection=False, subtract=True, selectregion=True,
                 **kwargs):
        """
        Fit and remove a polynomial from the spectrum.
        It will be saved in the variable "self.basespec"
        and the fit parameters will be saved in "self.order"

        Parameters
        ----------
        order: int
            Order of the polynomial to fit
        excludefit: bool
            If there is a spectroscopic line fit, you can automatically exclude
            the region with signal above some tolerance set by `exclusionlevel`
            (it works for absorption lines by using the absolute value of the
            signal)
        exclusionlevel: float
            The minimum value of the spectroscopic fit to exclude when fitting
            the baseline
        save: bool
            Write the baseline fit coefficients into the spectrum's header in the
            keywords BLCOEFnn
        interactive: bool
            Specify the include/exclude regions through the interactive plot
            window
        fit_original: bool
            Fit the original spectrum instead of the baseline-subtracted spectrum.
            If disabled, will overwrite the original data with the
            baseline-subtracted version.

            .. warning:: If this is set False, behavior of `unsubtract` may be unexpected
        fit_plotted_area: bool
            Will respect user-specified zoom (using the pan/zoom buttons)
            unless xmin/xmax have been set manually
        reset_selection: bool
            Reset the selected region to those specified by this command only
            (will override previous xmin/xmax settings)
        select_region: bool
            Run the region selection procedure?  If false, will leave
            'includemask' untouched

        baseline_fit_color: color name (string)
            [plot parameter]
            Color to plot the baseline
        clear_all_connections: bool
            [plot parameter]
            Disable any previous interactive sessions
        highlight_fitregion: bool
            [plot parameter]
            Highlight the selected region for baseline fitting (default green)

        """
        specfit = self.Spectrum.specfit
        self.order = order
        self.set_spectofit(fit_original=fit_original)
        self.OKmask = (self.spectofit==self.spectofit)
        if (isinstance(exclude, str) and (exclude == 'interactive')) or interactive:
            if np.all(self.includemask):
                specwarnings.warn("Include mask was uniformly True.  \n"
                                  "This has the effect of making the 'baseline'"
                                  " command do nothing.\n"
                                  "If this was not your intent, try"
                                  " spectrum.baseline(interactive=True,"
                                  "reset_selection=True)")
            self.start_interactive(clear_all_connections=clear_all_connections,
                                   debug=debug,
                                   reset_selection=reset_selection, **kwargs)
        else:
            if reset_selection:
                # this could be accomplished by passing reset=True to
                # selectregion, but that may cause
                # clashes [unclear; note added 2/12/2014]
                self.includemask[:] = True
                if exclude is not None:
                    log.warning("`reset_selection` was set, so `exclude` is being ignored.")
            elif selectregion:
                # must select region (i.e., exclude edges) AFTER setting 'positive'
                # include region
                # also, DON'T highlight here because it will be cleared
                self.selectregion(fit_plotted_area=fit_plotted_area,
                                  exclude=exclude, debug=debug, **kwargs)

            # have to do excludefit after selection, otherwise
            # fit_plotted_area=True overrides
            if excludefit and specfit.modelpars is not None:
                #vlo = self.Spectrum.plotter.specfit.modelpars[1] - 2*self.Spectrum.plotter.specfit.modelpars[2]
                #vhi = self.Spectrum.plotter.specfit.modelpars[1] + 2*self.Spectrum.plotter.specfit.modelpars[2]
                #exclude = [np.argmin(abs(self.Spectrum.xarr-vlo)),argmin(abs(self.Spectrum.xarr-vhi))]
                full_model = specfit.get_full_model(add_baseline=False)

                # only set additional FALSE
                self.includemask *= (abs(full_model) <
                                     exclusionlevel*np.abs(full_model).max())
                #import pdb; pdb.set_trace()

            log.debug("kwargs to button2 action = {0}".format(kwargs))


            self.button2action(fit_original=fit_original,
                               baseline_fit_color=baseline_fit_color,
                               debug=debug,
                               subtract=subtract,
                               LoudDebug=LoudDebug,
                               **kwargs)

            if highlight_fitregion:
                self.highlight_fitregion()
        if save:
            self.savefit()

    def set_spectofit(self, fit_original=True, fit_residuals=False):
        """
        Reset the spectrum-to-fit from the data
        """
        if fit_residuals:
            self.spectofit = self.Spectrum.specfit.residuals
        elif self.subtracted and fit_original: # add back in the old baseline
            self.spectofit = self.Spectrum.data+self.basespec
        else:
            self.spectofit = self.Spectrum.data.copy()

    def fit(self, powerlaw=None, order=None, includemask=None, spline=False,
            spline_sampling=10, spline_downsampler=np.median,
            xarr_fit_unit='pixels', **kwargs):
        """
        Run the fit and set `self.basespec`
        """
        log.debug("In baseline.fit, xarr_fit_unit={0}"
                  .format(xarr_fit_unit))

        self.powerlaw = powerlaw or self.powerlaw
        self.spline = spline
        self.order = order or self.order
        if includemask is not None:
            self.includemask = includemask

        spline_sampling = (1 if spline_sampling is None else spline_sampling)

        xarr = self.Spectrum.xarr
        if spline:
            self.basespec = _spline(self.spectofit,
                                    xarr=xarr,
                                    masktofit=self.includemask,
                                    sampling=spline_sampling,
                                    downsampler=spline_downsampler,
                                    order=self.order)
        else:
            # use a,b to keep line under 80 chars
            a,b = self._baseline(self.spectofit, xarr=xarr,
                                 err=self.Spectrum.error, order=self.order,
                                 masktoexclude=~self.includemask,
                                 powerlaw=self.powerlaw,
                                 xarr_fit_unit=xarr_fit_unit,
                                 **kwargs)
            self.basespec, self.baselinepars = (a,b)

            self.set_basespec_frompars(xarr_fit_unit=xarr_fit_unit)

    def button2action(self, event=None, debug=False, subtract=True,
                      powerlaw=None, fit_original=False,
                      spline=False,
                      spline_sampling=None,
                      spline_downsampler=np.median,
                      baseline_fit_color='orange', **kwargs):
        """
        Do the baseline fitting and save and plot the results.

        Note that powerlaw fitting will only consider positive data.
        """
        if debug:
            print("Button 2/3 Baseline.  Subtract={0}".format(subtract))
        if self.subtracted:
            self.unsubtract()

        # interactive guesspeakwidth passes this on; it should be excised
        if 'nwidths' in kwargs:
            kwargs.pop('nwidths')

        if powerlaw is None:
            powerlaw = self.powerlaw

        self.clear_highlights()

        self._xfit_units = self.Spectrum.xarr.unit

        log.debug("Fitting baseline: powerlaw={0} "
                  "spline={1}".format(powerlaw, spline))
        self.fit(powerlaw=powerlaw, includemask=self.includemask,
                 order=self.order, spline=spline,
                 spline_sampling=spline_sampling,
                 spline_downsampler=spline_downsampler,
                 **kwargs)

        if subtract:
            if self.subtracted and fit_original:
                # use the spectrum with the old baseline added in (that's what we fit to)
                self.Spectrum.data = self.spectofit - self.basespec
            else:
                self.Spectrum.data -= self.basespec
            self.subtracted = True
        else:
            if self.subtracted:
                self.unsubtract()
            self.subtracted = False

        if self.Spectrum.plotter.axis is not None:
            log.debug("Plotting baseline")
            if event is not None:
                # preserve frame if fitting interactively
                kwargs.update({'use_window_limits':True})
            self.plot_baseline(baseline_fit_color=baseline_fit_color, **kwargs)

        # disconnect interactive window (and more importantly, reconnect to
        # original interactive cmds)
        self.clear_all_connections()

        if hasattr(self.Spectrum,'header'):
            history.write_history(self.Spectrum.header,
                                  "BASELINE order=%i pars=%s"
                                  % (self.order, ",".join([str(s) for s in
                                                           self.baselinepars]))
                                  + "(powerlaw)" if self.powerlaw else "")

    def set_basespec_frompars(self, baselinepars=None, xarr_fit_unit=None):
        """
        Set the baseline spectrum based on the fitted parameters

        Parameters
        ----------
        baselinepars : list
            Optional list of fit parameters, e.g. a list of polynomial
            coefficients
        xarr_fit_unit : None or 'pixels' or 'native' or unit
            The units that were used in the baseline fit
        """
        if xarr_fit_unit == 'pixels':
            xarrconv = np.arange(self.Spectrum.data.size)
        elif xarr_fit_unit not in (None, 'native'):
            xarrconv = self.Spectrum.xarr.as_unit(xarr_fit_unit)
        else:
            xarrconv = self.Spectrum.xarr

        log.debug("xarr_fit_unit: {0}".format(xarr_fit_unit))

        self.basespec = self.get_model(xarr=xarrconv,
                                       baselinepars=baselinepars)

    def get_model(self, xarr=None, baselinepars=None):
        # create the full baseline spectrum...
        fit_units = self._xfit_units
        if xarr is None:
            xarr = self.Spectrum.xarr.as_unit(fit_units)
        if baselinepars is None:
            baselinepars = self.baselinepars
        if baselinepars is None: # still...
            return 0 # no baseline has been computed
        if self.powerlaw:
            result = (baselinepars[0]*(xarr)**(-baselinepars[1])).squeeze()
            if np.any(~np.isfinite(result)):
                raise ValueError("Infinite value encountered in power-law model")
            return result
        else:
            return np.poly1d(baselinepars)(xarr)


    def button3action(self, *args, **kwargs):
        """
        Wrapper - same as button2action, but with subtract=False
        """
        if 'subtract' in kwargs:
            kwargs.pop('subtract')

        return self.button2action(*args, subtract=False, **kwargs)

    def plot_baseline(self, annotate=True, baseline_fit_color='orange',
                      use_window_limits=None, linewidth=1, alpha=0.75,
                      plotkwargs={}, **kwargs):
        """
        Overplot the baseline fit

        Parameters
        ----------
        annotate : bool
            Display the fit parameters for the best-fit baseline on the
            top-left of the plot
        baseline_fit_color : matplotlib color
            What color to use for overplotting the line
            (default is slightly transparent orange)
        use_window_limits : None or bool
            Keep the current window or expand the plot limits?  If left as None,
            will use `self.use_window_limits`

        Other Parameters
        ----------------
        linewidth : number
        alpha : float [0-1]
        plotkwargs : dict
            Are passed to matplotlib's plot function
        """
        log.debug("At the beginning of plotbaseline, xmin, xmax = {0},{1}"
                  " and viewlim={2}"
                  .format(self.Spectrum.plotter.xmin,
                          self.Spectrum.plotter.xmax,
                          self.Spectrum.plotter.axis.viewLim))

        # clear out the errorplot.  This should not be relevant...
        if self.Spectrum.plotter.errorplot is not None:
            for p in self.Spectrum.plotter.errorplot:
                if isinstance(p,matplotlib.collections.PolyCollection):
                    if p in self.Spectrum.plotter.axis.collections:
                        try:
                            self.Spectrum.plotter.axis.lines.remove(p)
                        except AttributeError:
                            try:
                                p.remove()
                            except Exception as ex:
                                pass
                if isinstance(p,matplotlib.lines.Line2D):
                    if p in self.Spectrum.plotter.axis.lines:
                        try:
                            self.Spectrum.plotter.axis.lines.remove(p)
                        except AttributeError:
                            try:
                                p.remove()
                            except Exception as ex:
                                pass

        # if we subtract the baseline, replot the now-subtracted data with rescaled Y axes
        if self.subtracted:
            if self.Spectrum.plotter.axis is not None:
                for p in self.Spectrum.plotter.axis.lines:
                    try:
                        self.Spectrum.plotter.axis.lines.remove(p)
                    except AttributeError:
                        try:
                            p.remove()
                        except Exception as ex:
                            pass
            plotmask = self.OKmask*False # include nothing...
            plotmask[self.xmin:self.xmax] = self.OKmask[self.xmin:self.xmax] # then include everything OK in range
            self.Spectrum.plotter.ymin = abs(self.Spectrum.data[plotmask].min())*1.1*np.sign(self.Spectrum.data[plotmask].min())
            self.Spectrum.plotter.ymax = abs(self.Spectrum.data[plotmask].max())*1.1*np.sign(self.Spectrum.data[plotmask].max())
            # don't change the zoom (by default)!
            uwl = use_window_limits if use_window_limits is not None else self.use_window_limits
            self.Spectrum.plotter.plot(use_window_limits=uwl)
        else: # otherwise just overplot the fit
            self.Spectrum.plotter.axis.set_autoscale_on(False)
            for p in self._plots:
                # remove the old baseline plots
                if p in self.Spectrum.plotter.axis.lines:
                    try:
                        self.Spectrum.plotter.axis.lines.remove(p)
                    except AttributeError:
                        try:
                            p.remove()
                        except Exception as ex:
                            pass
            self._plots += self.Spectrum.plotter.axis.plot(
                    self.Spectrum.xarr,
                    self.basespec,
                    color=baseline_fit_color,
                    linewidth=linewidth,
                    alpha=alpha,
                    **plotkwargs)

        if annotate:
            self.annotate() # refreshes automatically
        elif self.Spectrum.plotter.autorefresh:
            self.Spectrum.plotter.refresh()

    def unsubtract(self, replot=True, preserve_limits=True):
        """
        Restore the spectrum to "pristine" state (un-subtract the baseline)

        *replot* [ True ]
            Re-plot the spectrum?  (only happens if unsubtraction proceeds,
            i.e. if there was a baseline to unsubtract)

        *preserve_limits* [ True ]
            Preserve the current x,y limits
        """
        if self.subtracted:
            self.Spectrum.data += self.basespec
            self.subtracted = False

            if replot and self.Spectrum.plotter.axis is not None:
                kwargs = self.Spectrum.plotter.plotkwargs
                kwargs.update({'use_window_limits':preserve_limits})
                if preserve_limits:
                    # need to manually stash the window limits because the plot
                    # is going to be cleared, and we need it to be cleared
                    self.Spectrum.plotter._stash_window_limits()
                    window_limits = self.Spectrum.plotter._window_limits
                self.Spectrum.plotter(**kwargs)

                if preserve_limits:
                    self.Spectrum.plotter._window_limits = window_limits
                    self.Spectrum.plotter._reset_to_stashed_limits()
        else:
            print("Baseline wasn't subtracted; not unsubtracting.")

    def annotate(self, loc='upper left', fontsize=10):
        if self.powerlaw:
            bltext = "bl: $y=%6.3g\\times(x)^{-%6.3g}$" % (self.baselinepars[0],
                                                           self.baselinepars[1])
        elif self.spline:
            bltext = "bl: spline"
        else:
            bltext = "bl: $y=$"+"".join(["$%+6.3gx^{%i}$" % (f,self.order-i)
                                         for i,f in
                                         enumerate(self.baselinepars)])
        #self.blleg = text(xloc,yloc     ,bltext,transform = self.Spectrum.plotter.axis.transAxes)
        self.clearlegend()
        pl = matplotlib.collections.CircleCollection([0],edgecolors=['k'])
        self.blleg = self.Spectrum.plotter.axis.legend((pl,), (bltext,),
                                                       loc=loc,
                                                       markerscale=0.001,
                                                       borderpad=0.1,
                                                       handlelength=0.1,
                                                       handletextpad=0.1,
                                                       frameon=False,
                                                       fontsize=fontsize)
        self.Spectrum.plotter.axis.add_artist(self.blleg)
        if self.Spectrum.plotter.autorefresh:
            self.Spectrum.plotter.refresh()

    def clearlegend(self):
        if self.blleg is not None:
            self.blleg.set_visible(False)
            if self.blleg in self.Spectrum.plotter.axis.artists:
                if hasattr(self.Spectrum.plotter.axis.artists, 'remove'):
                    self.Spectrum.plotter.axis.artists.remove(self.blleg)
                else:
                    self.blleg.remove()
        if self.Spectrum.plotter.autorefresh: self.Spectrum.plotter.refresh()

    def savefit(self):
        if self.baselinepars is not None and hasattr(self.Spectrum,'header'):
            for ii,p in enumerate(self.baselinepars):
                self.Spectrum.header['BLCOEF%0.2i' % (ii)] = (p,"Baseline power-law best-fit coefficient x^%i" % (self.order-ii-1))

    def _baseline(self, spectrum, xarr=None, err=None,
                  order=1, quiet=True, masktoexclude=None, powerlaw=False,
                  xarr_fit_unit='pixels', LoudDebug=False, renormalize='auto',
                  pguess=None, zeroerr_is_OK=True, spline=False, **kwargs):
        """
        Fit a baseline/continuum to a spectrum

        Parameters
        ----------
        masktoexclude : boolean array
            True: will not be fit
            False: will be fit
            *if ALL are True, ALL pixels will be fit - just with silly weights*
        """

        #if xmin == 'default':
        #    if order <= 1 and mask is None: xmin = np.floor( spectrum.shape[-1]*0.1 )
        #    else: xmin = 0
        #elif xmin is None:
        #    xmin = 0
        #if xmax == 'default':
        #    if order <= 1 and mask is None: xmax = np.ceil( spectrum.shape[-1]*0.9 )
        #    else: xmax = spectrum.shape[-1]
        #elif xmax is None:
        #    xmax = spectrum.shape[-1]

        if xarr is None or xarr_fit_unit == 'pixels':
            if powerlaw:
                # start at 1, not 0, for power-law in particular
                xarrconv = np.arange(1, spectrum.size+1)
            else:
                xarrconv = np.arange(spectrum.size)
        elif xarr_fit_unit not in (None, 'native'):
            xarrconv = xarr.as_unit(xarr_fit_unit).value
        else:
            xarrconv = xarr
        log.debug("xarrconv = {0}".format(xarrconv))
        log.debug("xarr_fit_unit = {0}".format(xarr_fit_unit))


        # A good alternate implementation of masking is to only pass mpfit the data
        # that is unmasked.  That would require some manipulation above...
        if err is None:
            err = np.ma.ones(spectrum.shape)
        else:
            # don't overwrite error
            err = np.ma.copy(err)
            # assume anything with 0 error is GOOD
            if zeroerr_is_OK:
                err[err == 0] = 1.
            else: # flag it out!
                err[err == 0] = 1e10
                err.mask[err == 0] = True


        #err[:xmin] = 1e10
        #err[xmax:] = 1e10
        if masktoexclude is not None:
            if masktoexclude.dtype.name != 'bool':
                masktoexclude = masktoexclude.astype('bool')
            err[masktoexclude] = 1e10
            err.mask |= masktoexclude
            if LoudDebug:
                print("In _baseline: %i points masked out" % masktoexclude.sum())
        if (spectrum!=spectrum).sum() > 0:
            print("There is an error in baseline: some values are NaN")
            import pdb; pdb.set_trace()

        OK = (~masktoexclude) & (~err.mask)
        if powerlaw:
            # for powerlaw fitting, only consider positive data
            OK &= spectrum > 0
            if pguess is None:
                pguess = [np.median(spectrum[OK]) * np.median(xarrconv[OK]), 1.0]
            if LoudDebug:
                print("_baseline powerlaw Guesses: ",pguess)

            parinfo = [{}, {'limited':(True, False), 'limits':(0,0)}]

            if any(xarrconv == 0) or any(~np.isfinite(xarrconv)):
                raise ValueError("Invalid value in x-axis being used for fitting.")

            def mpfitfun(data, err):
                #def f(p,fjac=None): return [0,np.ravel(((p[0] * (xarrconv[OK]-p[2])**(-p[1]))-data)/err)]
                # Logarithmic fitting:
                def f(params, fjac=None):
                    for pp in params:
                        if not np.isfinite(pp):
                            raise ValueError("Non-finite parameter value passed to baseine parameters")
                    #return [0,
                    #        np.ravel( (np.log10(data) - np.log10(p[0]) + p[1]*np.log10(xarrconv[OK]/p[2])) / (err/data) )
                    #        ]
                    rslt = np.ravel((data -
                                     self.get_model(xarr=xarrconv[OK],
                                                    baselinepars=params))
                                    / (err/data))
                    if np.any(np.isnan(rslt)):
                        raise ValueError("NaN in result")
                    return [0, rslt]
                return f
        else:
            pguess = [0]*(order+1)
            if LoudDebug:
                print("_baseline Guesses: ",pguess)

            parinfo = None

            def mpfitfun(data,err):
                def f(p,fjac=None):
                    return [0,np.ravel((np.poly1d(p)(xarrconv[OK])-data)/err)]
                return f
        #scalefactor = 1.0
        #if renormalize in ('auto',True):
        #    datarange = spectrum.max() - spectrum.min()
        #    if abs(datarange) < 1e-9 or abs(datarange) > 1e9:
        #        scalefactor = np.median(np.abs(self.spectrum))
        #        print "BASELINE: Renormalizing data by factor %e to improve fitting procedure" % scalefactor
        #        spectrum /= scalefactor
        #        err /= scalefactor

        # error checking:
        shapeOK = (xarrconv.shape == spectrum.shape and xarrconv.shape ==
                   err.shape and OK.shape == xarrconv.shape)
        if not shapeOK:
            raise ValueError("The shape of the fitted spectra are not"
                             " identical.  This is a major error, please "
                             "report it as an Issue: "
                             "https://github.com/pyspeckit/pyspeckit/issues")

        mp = mpfit.mpfit(mpfitfun(np.array(spectrum[OK]), np.array(err[OK])), xall=pguess,
                         parinfo=parinfo,
                         quiet=quiet)
        if np.isnan(mp.fnorm):
            raise ValueError("chi^2 is NAN in baseline fitting")
        self.mpfit_status = models.mpfit_messages[mp.status]
        log.debug("mpfit message: {0}".format(self.mpfit_status))
        fitp = mp.params
        if powerlaw:
            bestfit = (fitp[0]*(xarrconv)**(-fitp[1])).squeeze()
        else:
            bestfit = np.poly1d(fitp)(xarrconv).squeeze()

        return bestfit,fitp


    def crop(self,x1pix,x2pix):
        """
        When spectrum.crop is called, this must be too
        """
        self.basespec = self.basespec[x1pix:x2pix]
        self.includemask = self.includemask[x1pix:x2pix]
        self.OKmask = self.OKmask[x1pix:x2pix]

    def downsample(self,factor):
        self.basespec = self.basespec[::factor]
        self.includemask = self.includemask[::factor]
        self.OKmask = self.OKmask[::factor]

    def copy(self, parent=None):
        """
        Create a copy of the baseline fit

        [ parent ]
            A spectroscopic axis instance that is the parent of the specfit
            instance.  This needs to be specified at some point, but defaults
            to None to prevent overwriting a previous plot.
        """

        newbaseline = copy.copy(self)
        newbaseline.Spectrum = parent
        newbaseline.OKmask = copy.copy( self.OKmask )
        newbaseline.basespec = copy.copy( self.basespec )
        newbaseline.baselinepars = copy.copy( self.baselinepars )
        newbaseline.includemask = self.includemask.copy()

        if parent is not None:
            newbaseline.Spectrum.plotter = parent.plotter
        else:
            newbaseline.Spectrum.plotter = None

        return newbaseline

def _spline(data, xarr=None, masktofit=None, order=3, sampling=10,
            downsampler=np.median, append_endpoints=True):
    from scipy.interpolate import UnivariateSpline

    assert masktofit.shape == data.shape

    if masktofit is None:
        masktofit = np.isfinite(data)
        if not any(masktofit):
            log.warning("All data was infinite or NaN")

    if xarr is None:
        xarr = np.arange(data.size, dtype=data.dtype)

    if downsampler is not None:
        ngood = np.count_nonzero(masktofit)
        if ngood == 0:
            log.warning("Fitting all data with spline")
            masktofit[:] = True
            ngood = masktofit.size
        endpoint = ngood - (ngood % sampling)
        yfit = downsampler([data[masktofit][ii:endpoint:sampling]
                            for ii in range(sampling)], axis=0)
        xfit = xarr[masktofit][sampling/2:endpoint:sampling]
        if append_endpoints:
            xfit = np.array([xarr[0]-(xarr[1]-xarr[0])] +
                            xfit.tolist() +
                            [xarr[-1]+(xarr[1]-xarr[0])])
            yfit = np.array([yfit[0]] + yfit.tolist() + [yfit[-1]])
        inds = np.argsort(xfit)
        spl = UnivariateSpline(xfit[inds],
                               yfit[inds],
                               k=order,
                               s=0)
    else:
        xfit = xarr[masktofit]
        yfit = data[masktofit]
        inds = np.argsort(xfit)
        spl = UnivariateSpline(xfit[inds],
                               yfit[inds],
                               k=order,
                               s=sampling)

    baseline_fitted = spl(xarr)
    if np.any(np.isnan(baseline_fitted)):
        log.error("NaNs in baseline.")
        import ipdb; ipdb.set_trace()
    return baseline_fitted

