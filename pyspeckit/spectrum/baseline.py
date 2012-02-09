import numpy as np
import matplotlib
from ..config import mycfg
from ..config import ConfigDescriptor as cfgdec
import interactive

interactive_help_message = """
Left-click twice to select or add to the baseline fitting range.  Middle or
right click to disconnect and perform the fit.
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
        self.specplotter = Spectrum.plotter
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
        self._plots = []

    @cfgdec
    def __call__(self, order=1, excludefit=False, save=True, exclude=None,
            exclusionlevel=0.01, interactive=False, debug=False,
            LoudDebug=False, fit_original=True, baseline_fit_color='orange',
            clear_all_connections=True, fit_plotted_area=True, highlight=False,
            reset_selection=False,
            **kwargs):
        """
        Fit and remove a polynomial from the spectrum.  
        It will be saved in the variable "self.basespec"
        and the fit parameters will be saved in "self.order"

        function baseline(spectrum,xarr=None,xmin=None,xmax=None,order=1,quiet=True,exclude=None):
        Subtract a baseline from a spectrum
        If xmin,xmax are not specified, defaults to ignoring first and last 10% of spectrum

        exclude is a set of start/end indices to ignore when baseline fitting
        (ignored by setting error to infinite in fitting procedure)

        excludefit creates a mask based on the fitted gaussian model (assuming
        that it has a zero-height) using an exclusion level of (exclusionlevel)
        * the smallest gaussian peak that was fit

        fit_plotted_area means that a fit will only be attempted on the region
        currently shown in the plot.  Overridden by interactive=True or
        exclude='interactive'

        if fit_original is set, "basespec" is added back to the spectrum before
        fitting so you can run this procedure multiple times without losing
        information
        """
        specfit = self.Spectrum.specfit
        self.order = order
        fitp = np.zeros(self.order+1)
        if self.subtracted and fit_original: # add back in the old baseline
            self.spectofit = self.Spectrum.data+self.basespec
        else:
            self.spectofit = self.Spectrum.data.copy()
        self.OKmask = (self.spectofit==self.spectofit)
        if exclude == 'interactive' or interactive:
            self.start_interactive(clear_all_connections=clear_all_connections,
                    debug=debug, **kwargs)
        else:
            if excludefit and specfit.modelpars is not None:
                #vlo = self.specplotter.specfit.modelpars[1] - 2*self.specplotter.specfit.modelpars[2]
                #vhi = self.specplotter.specfit.modelpars[1] + 2*self.specplotter.specfit.modelpars[2]
                #exclude = [np.argmin(abs(self.Spectrum.xarr-vlo)),argmin(abs(self.Spectrum.xarr-vhi))]
                specfit.fullsizemodel() # make sure the spectrum is the right size
                if reset_selection:
                    self.includemask = abs(specfit.model) < exclusionlevel*min(abs(np.array(specfit.modelpars[0::3])))
                else:
                    # only set additional FALSE
                    self.includemask *= abs(specfit.model) < exclusionlevel*min(abs(np.array(specfit.modelpars[0::3])))
            elif reset_selection:
                self.includemask[:] = True
            # must select region (i.e., exclude edges) AFTER setting 'positive' include region
            # also, DON'T highlight here because it will be cleared
            self.selectregion(fit_plotted_area=fit_plotted_area, **kwargs)
            self.button3action(exclude=exclude, 
                    fit_original=fit_original,
                    baseline_fit_color=baseline_fit_color, **kwargs)
            if highlight: self.highlight_fitregion()
        if save: self.savefit()

    def button3action(self, event=None, debug=False, subtract=True,
            fit_original=False, powerlaw=False, baseline_fit_color='orange',
            exclude=None, **kwargs):
        """
        Do the baseline fitting and save and plot the results.

        Can specify a region to exclude using velocity units or pixel units
        """
        if debug: print "Button 3 Baseline"
        if self.subtracted:
            self.unsubtract()

        self.clear_highlights()

        xarr_fit_units = self.Spectrum.xarr.units

        # Exclude keyword-specified excludes.  Assumes exclusion in current X array units
        if exclude is not None and len(exclude) % 2 == 0:
            for x1,x2 in zip(exclude[::2],exclude[1::2]):
                x1 = self.Spectrum.xarr.x_to_pix(x1)
                x2 = self.Spectrum.xarr.x_to_pix(x2)
                self.includemask[x1:x2] = False

        self.basespec, self.baselinepars = self._baseline(
                self.spectofit,
                xarr=self.Spectrum.xarr,
                err=self.Spectrum.error,
                order=self.order, 
                mask=(True-self.includemask),
                powerlaw=powerlaw,
                xarr_fit_units=xarr_fit_units,
                **kwargs)
        # create the full baseline spectrum...
        if powerlaw:
            self.powerlaw = True
            #self.basespec = (self.baselinepars[0]*(self.Spectrum.xarr.as_unit(xarr_fit_units)-self.baselinepars[2])**(-self.baselinepars[1])).squeeze()
            self.basespec = (self.baselinepars[0]*(self.Spectrum.xarr.as_unit(xarr_fit_units)/self.baselinepars[2])**(-self.baselinepars[1])).squeeze()
        else:
            self.powerlaw = False
            self.basespec = np.poly1d(self.baselinepars)(self.Spectrum.xarr.as_unit(xarr_fit_units))
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

        if self.specplotter.axis is not None:
            self.plot_baseline(baseline_fit_color=baseline_fit_color, **kwargs)

        # disconnect interactive window
        self.clear_all_connections()

    button2action = button3action

    def plot_baseline(self, annotate=True, baseline_fit_color='orange',
            **kwargs):
        """
        Overplot the baseline fit
        """

        # clear out the errorplot.  This should not be relevant...
        if self.specplotter.errorplot is not None: 
            for p in self.specplotter.errorplot:
                if isinstance(p,matplotlib.collections.PolyCollection):
                    if p in self.specplotter.axis.collections: 
                        self.specplotter.axis.collections.remove(p)
                if isinstance(p,matplotlib.lines.Line2D):
                    if p in self.specplotter.axis.lines: 
                        self.specplotter.axis.lines.remove(p)

        # if we subtract the baseline, replot the now-subtracted data with rescaled Y axes
        if self.subtracted:
            if self.specplotter.axis is not None: 
                for p in self.specplotter.axis.lines:
                    self.specplotter.axis.lines.remove(p)
            plotmask = self.OKmask*False # include nothing...
            plotmask[self.xmin:self.xmax] = self.OKmask[self.xmin:self.xmax] # then include everything OK in range
            self.specplotter.ymin = abs(self.Spectrum.data[plotmask].min())*1.1*np.sign(self.Spectrum.data[plotmask].min())
            self.specplotter.ymax = abs(self.Spectrum.data[plotmask].max())*1.1*np.sign(self.Spectrum.data[plotmask].max())
            self.specplotter.plot()
        else: # otherwise just overplot the fit
            self.specplotter.axis.set_autoscale_on(False)
            for p in self._plots:
                # remove the old baseline plots
                if p in self.specplotter.axis.lines:
                    self.specplotter.axis.lines.remove(p)
            self._plots += self.specplotter.axis.plot(self.Spectrum.xarr,self.basespec,color=baseline_fit_color)

        if annotate: self.annotate() # refreshes automatically
        elif self.specplotter.autorefresh: self.specplotter.refresh()

    def unsubtract(self):
        if self.subtracted:
            self.Spectrum.data += self.basespec
            self.subtracted = False
        else: 
            print "Baseline wasn't subtracted; not unsubtracting."

    def annotate(self,loc='upper left'):
        if self.powerlaw:
            #bltext = "bl: $y=%6.3g\\times(x-%6.3g)^{-%6.3g}$" % (self.baselinepars[0],self.baselinepars[2],self.baselinepars[1])
            bltext = "bl: $y=%6.3g\\times(x/%6.3g)^{-%6.3g}$" % (self.baselinepars[0],self.baselinepars[2],self.baselinepars[1])
        else:
            bltext = "bl: $y=$"+"".join(["$%+6.3gx^{%i}$" % (f,self.order-i)
                for i,f in enumerate(self.baselinepars)])
        #self.blleg = text(xloc,yloc     ,bltext,transform = self.specplotter.axis.transAxes)
        self.clearlegend()
        pl = matplotlib.collections.CircleCollection([0],edgecolors=['k'])
        self.blleg = self.specplotter.axis.legend(
                (pl,),
                (bltext,),loc=loc,markerscale=0.001,
                borderpad=0.1, handlelength=0.1, handletextpad=0.1, frameon = False
                )
        self.specplotter.axis.add_artist(self.blleg)
        if self.specplotter.autorefresh: self.specplotter.refresh()
  
    def clearlegend(self):
        if self.blleg is not None: 
            self.blleg.set_visible(False)
            if self.blleg in self.specplotter.axis.artists:
                self.specplotter.axis.artists.remove(self.blleg)
        if self.specplotter.autorefresh: self.specplotter.refresh()

    def savefit(self):
        if self.baselinepars is not None and hasattr(self.Spectrum,'header'):
            for ii,p in enumerate(self.baselinepars):
                self.Spectrum.header.update('BLCOEF%0.2i' % (ii),p,comment="Baseline power-law best-fit coefficient x^%i" % (self.order-ii-1))

    def _baseline(self, spectrum, xarr=None, err=None, xmin='default', xmax='default',
            order=1, quiet=True, mask=None, powerlaw=False,
            xarr_fit_units='pixels', LoudDebug=False, renormalize='auto',
            zeroerr_is_OK=True, **kwargs):
        """
        Subtract a baseline from a spectrum
        If xmin,xmax are not specified, defaults to ignoring first and last 10% of spectrum
        *unless* order > 1, in which case ignoring the ends tends to cause strange effects

        EDIT Nov 21, 2011: DO NOT exclude ends by chopping them!  This results
        in bad fits when indexing by pixels
        """

        if xmin == 'default':
            if order <= 1 and mask is None: xmin = np.floor( spectrum.shape[-1]*0.1 )
            else: xmin = 0
        elif xmin is None:
            xmin = 0
        if xmax == 'default':
            if order <= 1 and mask is None: xmax = np.ceil( spectrum.shape[-1]*0.9 )
            else: xmax = spectrum.shape[-1]
        elif xmax is None:
            xmax = spectrum.shape[-1]
        
        if xarr is None:
            xarr = np.indices(spectrum.shape).squeeze()


        # A good alternate implementation of masking is to only pass mpfit the data
        # that is unmasked.  That would require some manipulation above...
        if err is None:
            err = np.ones(spectrum.shape)
        else:
            # don't overwrite error
            err = err.copy()
            # assume anything with 0 error is GOOD
            if zeroerr_is_OK:
                err[err == 0] = 1.
            else: # flag it out!
                err[err == 0] = 1e10


        err[:xmin] = 1e10
        err[xmax:] = 1e10
        if mask is not None:
            if mask.dtype.name != 'bool': mask = mask.astype('bool')
            err[mask] = 1e10
            if LoudDebug: print "In _baseline: %i points masked out" % mask.sum()
        if (spectrum!=spectrum).sum() > 0:
            print "There is an error in baseline: some values are NaN"
            import pdb; pdb.set_trace()

        #xarrconv = xarr[xmin:xmax].as_unit(xarr_fit_units)
        OK = True-mask
        OK[:xmin] = False
        OK[xmax:] = False
        xarrconv = xarr.as_unit(xarr_fit_units)
        if powerlaw:
            pguess = [np.median(spectrum[OK]),2,xarrconv[OK][0]-1]
            if LoudDebug: print "_baseline powerlaw Guesses: ",pguess

            def mpfitfun(data,err):
                #def f(p,fjac=None): return [0,np.ravel(((p[0] * (xarrconv[OK]-p[2])**(-p[1]))-data)/err)]
                # Logarithmic fitting:
                def f(p,fjac=None):
                    return [0,
                            np.ravel( (np.log10(data) - np.log10(p[0]) + p[1]*np.log10(xarrconv[OK]/p[2])) / (err/data) )
                            ]
                return f
        else:
            pguess = [0]*(order+1)
            if LoudDebug: print "_baseline Guesses: ",pguess

            def mpfitfun(data,err):
                def f(p,fjac=None): return [0,np.ravel((np.poly1d(p)(xarrconv[OK])-data)/err)]
                return f
        #scalefactor = 1.0
        #if renormalize in ('auto',True):
        #    datarange = spectrum.max() - spectrum.min()
        #    if abs(datarange) < 1e-9 or abs(datarange) > 1e9:
        #        scalefactor = np.median(np.abs(self.spectrum))
        #        print "BASELINE: Renormalizing data by factor %e to improve fitting procedure" % scalefactor
        #        spectrum /= scalefactor
        #        err /= scalefactor

        if LoudDebug:
            print "In _baseline, xmin: %i, xmax: %i, len(spectrum)=%i, len(err)=%i" % (xmin,xmax,len(spectrum),len(err))

        import pyspeckit.mpfit as mpfit
        mp = mpfit.mpfit(mpfitfun(spectrum[OK],err[OK]),xall=pguess,quiet=quiet) # mpfit doesn't need to take kwargs, I think ,**kwargs)
        if np.isnan(mp.fnorm):
            raise ValueError("chi^2 is NAN in baseline fitting")
        fitp = mp.params
        if powerlaw:
            #bestfit = (fitp[0]*(xarrconv-fitp[2])**(-fitp[1])).squeeze()
            bestfit = (fitp[0]*(xarrconv/fitp[2])**(-fitp[1])).squeeze()
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

