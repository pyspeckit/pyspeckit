import matplotlib
import matplotlib.cbook as mpcb
import matplotlib.pyplot as pyplot
import numpy as np
from ..config import mycfg
from ..config import ConfigDescriptor as cfgdec
import units
import models
from pyspeckit.specwarnings import warn
import interactive
import inspect
import copy

class Registry(object):
    """
    This class is a simple wrapper to prevent fitter properties from being globals
    """

    def __init__(self):
        self.npars = {}
        self.multifitters = {}
        self.singlefitters = {}
        self.fitkeys = {}
        self.associatedkeys = {}

        self.interactive_help_message = """
        Left-click or hit 'p' twice to select (/p/ick) a fitting range.  You
        can e/x/clude or /r/emove parts of the spectrum by hitting 'x' or 'r'
        twice.  Then middle-click or hit 'm' twice to select (/m/ark) a peak
        and width.  When you're done, right-click or hit 'd' to perform the fit
        and disconnect the mouse and keyboard (/d/isconnect because you're
        /d/one).
        '?' will print this help message again.
        You can select different fitters to use with the interactive fitting routine.
        The default is gaussian ('g')
        """

    def add_fitter(self, name, function, npars, multisingle='single',
        override=False, key=None):
        ''' 
        Register a fitter function.

        Required Arguments:

            *name*: [ string ]
                The fit function name. 

            *function*: [ function ]
                The fitter function.  Single-fitters should take npars + 1 input
                parameters, where the +1 is for a 0th order baseline fit.  They
                should accept an X-axis and data and standard fitting-function
                inputs (see, e.g., gaussfitter).  Multi-fitters should take N *
                npars, but should also operate on X-axis and data arguments.

            *npars*: [ int ]
                How many parameters does the function being fit accept?

        Optional Keyword Arguments:

            *multisingle*: [ 'multi' | 'single' ] 
                Is the function a single-function fitter (with a background), or
                does it allow N copies of the fitting function?

            *override*: [ True | False ]
                Whether to override any existing type if already present.

            *key*: [ char ]
                Key to select the fitter in interactive mode
        '''


        if multisingle == 'single':
            if not name in self.singlefitters or override:
                self.singlefitters[name] = function
        elif multisingle == 'multi':
            if not name in self.multifitters or override:
                self.multifitters[name] = function
        elif name in self.singlefitters or name in self.multifitters:
            raise Exception("Fitting function %s is already defined" % name)

        if key is not None:
            self.fitkeys[key] = name
            self.interactive_help_message += "\n'%s' - select fitter %s" % (key,name)
        self.npars[name] = npars
        self.associated_keys = dict(zip(self.fitkeys.values(),self.fitkeys.keys()))

default_Registry = Registry()
default_Registry.add_fitter('ammonia',models.ammonia_model(multisingle='multi'),6,multisingle='multi',key='a')
default_Registry.add_fitter('ammonia_tau',models.ammonia_model_vtau(multisingle='multi'),6,multisingle='multi')
# not implemented default_Registry.add_fitter(Registry,'ammonia',models.ammonia_model(multisingle='single'),6,multisingle='single',key='A')
default_Registry.add_fitter('formaldehyde',models.formaldehyde_fitter,3,multisingle='multi',key='F') # CAN'T USE f!  reserved for fitting
default_Registry.add_fitter('formaldehyde',models.formaldehyde_vheight_fitter,3,multisingle='single')
default_Registry.add_fitter('gaussian',models.gaussian_fitter(multisingle='multi'),3,multisingle='multi',key='g')
default_Registry.add_fitter('gaussian',models.gaussian_fitter(multisingle='single'),3,multisingle='single')
default_Registry.add_fitter('voigt',models.voigt_fitter(multisingle='multi'),4,multisingle='multi',key='v')
default_Registry.add_fitter('voigt',models.voigt_fitter(multisingle='single'),4,multisingle='single')
default_Registry.add_fitter('lorentzian',models.lorentzian_fitter(multisingle='multi'),3,multisingle='multi',key='L')
default_Registry.add_fitter('lorentzian',models.lorentzian_fitter(multisingle='single'),3,multisingle='single')
default_Registry.add_fitter('hill5',models.hill5infall.hill5_fitter,5,multisingle='multi')
default_Registry.add_fitter('hcn',models.hcn.hcn_vtau_fitter,4,multisingle='multi')


class Specfit(interactive.Interactive):

    def __init__(self, Spectrum, Registry=None):
        super(Specfit, self).__init__(Spectrum, interactive_help_message=Registry.interactive_help_message)
        self.model = None
        self.modelpars = None
        self.modelerrs = None
        self.modelplot = []
        self.modelcomponents = None
        self._plotted_components = []
        self.npeaks = 0
        #self.nclicks_b1 = 0
        #self.nclicks_b2 = 0
        #self.xmin = 0
        #self.xmax = Spectrum.data.shape[0]
        self.button2action = self.guesspeakwidth
        self.guesses = []
        self.click = 0
        self.fitkwargs = {}
        self.auto = False
        #self.Spectrum = Spectrum
        self.Spectrum.plotter = self.Spectrum.plotter
        self.fitleg=None
        self.residuals=None
        self.setfitspec()
        self.fittype = 'gaussian'
        self.measurements = None
        self.vheight=False # vheight must be a boolean, can't be none
        self.Registry = Registry
        self.autoannotate = mycfg['autoannotate']
        #self.seterrspec()
        
    @cfgdec
    def __call__(self, interactive=False, usemoments=True,
            clear_all_connections=True, debug=False, multifit=False,
            guesses=None, save=True, fittype='gaussian', annotate=None,
            show_components=None,
            verbose=True, clear=True, vheight=None, **kwargs):
        """
        Fit gaussians (or other model functions) to a spectrum

        guesses = [height,amplitude,center,width]

        If you pass interactive=True, you can fit n gaussians using the mouse
        and/or keyboard:
            Left click or 'p': Set fitting region.  Two clicks sets the region,
                which will be highlighted.
            Middle click or 'm': Select peaks and peak widths.  The first click
                will mark a peak X-location and height with an X, the second
                click will mark the half-width-half-max location with a red
                line that represents the full-width-half-max
            Right click or 'd': Disconnect the plot and perform the fit.
        """

        if clear: self.clear()
        self.selectregion(verbose=verbose, debug=debug, **kwargs)
        for arg in ['xmin','xmax','xtype','reset']: 
            if arg in kwargs: kwargs.pop(arg)

        if self.fittype != fittype: self.fittype = fittype

        # multifit = True if the right number of guesses are passed
        if guesses is not None:
            if len(guesses) > 5:
                multifit = True

        self.npeaks = 0
        self.fitkwargs = kwargs
        if interactive:
            if self.Spectrum.plotter.axis is None:
                raise Exception("Interactive fitting requires a plotter.")
            # reset button count & guesses on every __call__
            self.nclicks_b1 = 0
            self.nclicks_b2 = 0
            self.guesses = []

            self.start_interactive(clear_all_connections=clear_all_connections, debug=debug, **kwargs)
        elif ((multifit or multifit is None) and self.fittype in self.Registry.multifitters) or guesses is not None:
            if guesses is None:
                print "You must input guesses when using multifit.  Also, baseline (continuum fit) first!"
                return
            else:
                self.guesses = guesses
                self.multifit(show_components=show_components, verbose=verbose, debug=debug)
        # SINGLEFITTERS SHOULD BE PHASED OUT
        elif self.fittype in self.Registry.singlefitters:
            #print "Non-interactive, 1D fit with automatic guessing"
            if (self.Spectrum.baseline.order is None and vheight is None) or vheight:
                self.Spectrum.baseline.order=0
                self.peakbgfit(usemoments=usemoments, show_components=show_components, annotate=annotate, debug=debug, **kwargs)
            else:
                self.peakbgfit(usemoments=usemoments, vheight=False,
                        height=0.0, annotate=annotate,
                        show_components=show_components, debug=debug, **kwargs)
            if self.Spectrum.plotter.autorefresh: self.Spectrum.plotter.refresh()
        else:
            if multifit:
                print "Can't fit with given fittype %s: it is not Registryed as a multifitter." % self.fittype
            else:
                print "Can't fit with given fittype %s: it is not Registryed as a singlefitter." % self.fittype
            return
        if save: self.savefit()


    def EQW(self, plot=False, plotcolor='g', annotate=False, alpha=0.5, loc='lower left'):
        """
        Returns the equivalent width (integral of "baseline" or "continuum"
        minus the spectrum) over the selected range
        """
        if np.median(self.Spectrum.baseline.basespec) == 0:
            raise ValueError("Baseline / continuum is zero: equivalent width is undefined.")
        elif np.median(self.Spectrum.baseline.basespec) < 0:
            if mycfg.WARN: warn( "WARNING: Baseline / continuum is negative: equivalent width is poorly defined." )
        diffspec = (self.Spectrum.baseline.basespec - self.Spectrum.data)
        dx = np.abs((self.Spectrum.xarr[self.xmax-1]-self.Spectrum.xarr[self.xmin]) / (self.xmax-self.xmin))
        sumofspec = diffspec[self.xmin:self.xmax].sum() * dx
        eqw = sumofspec / np.median(self.Spectrum.baseline.basespec)
        if plot:
            midpt_pixel = np.round((self.xmin+self.xmax)/2.0)
            midpt       = self.Spectrum.xarr[midpt_pixel]
            midpt_level = self.Spectrum.baseline.basespec[midpt_pixel]
            print "EQW plotting: ",midpt,midpt_pixel,midpt_level,eqw
            self.Spectrum.plotter.axis.fill_between(
                    [midpt-eqw/2.0,midpt+eqw/2.0],
                    [0,0],
                    [midpt_level,midpt_level],
                    color=plotcolor,
                    alpha=alpha,
                    label='EQW: %0.3g' % eqw)
            if annotate:
                self.Spectrum.plotter.axis.legend(
                        [(matplotlib.collections.CircleCollection([0],facecolors=[plotcolor],edgecolors=[plotcolor]))],
                        [('EQW: %0.3g' % eqw)], 
                        markerscale=0.01, borderpad=0.1, handlelength=0.1,
                        handletextpad=0.1, loc=loc)
            if self.Spectrum.plotter.autorefresh:
                self.Spectrum.plotter.refresh()
        return eqw

    def register_fitter(self,*args,**kwargs):
        """
        Register a model fitter
        """
        self.Registry.add_fitter(*args,**kwargs)

    register_fitter.__doc__ += Registry.add_fitter.__doc__
    
    def seterrspec(self,usestd=None,useresiduals=True):
        """
        Simple wrapper function to set the error spectrum; will either use the
        input spectrum or determine the error using the RMS of the residuals,
        depending on whether the residuals exist.
        """
        if self.residuals is not None and useresiduals: 
            self.errspec = np.ones(self.spectofit.shape[0]) * self.residuals.std()
        elif self.Spectrum.error is not None and not usestd:
            if (self.Spectrum.error == 0).all():
                if type(self.Spectrum.error) is np.ma.masked_array:
                    # force errspec to be a non-masked array of ones
                    self.errspec = self.Spectrum.error.data + 1
                else:
                    self.errspec = self.Spectrum.error + 1
            else:
                self.errspec = self.Spectrum.error
        else: self.errspec = np.ones(self.spectofit.shape[0]) * self.spectofit.std()

    def setfitspec(self):
        """
        Set the spectrum that will be fit.  This is primarily to remove NANs
        from consideration: if you simply remove the data from both the X-axis
        and the Y-axis, it will not be considered for the fit, and a linear
        X-axis is not needed for fitting.

        However, it may be possible to do this using masked arrays instead of
        setting errors to be 1e10....
        """
        self.spectofit = np.copy(self.Spectrum.data)
        if hasattr(self.Spectrum,'baseline'):
            if (self.Spectrum.baseline.subtracted is False 
                    and self.Spectrum.baseline.basespec is not None
                    and len(self.spectofit) == len(self.Spectrum.baseline.basespec)):
                self.spectofit -= self.Spectrum.baseline.basespec
        OKmask = (self.spectofit==self.spectofit)
        self.spectofit[(True-OKmask)] = 0
        self.seterrspec()
        self.errspec[(True-OKmask)] = 1e10
        if self.includemask is not None and (self.includemask.shape == self.errspec.shape):
            self.errspec[True - self.includemask] = 1e10

    def multifit(self, fittype=None, renormalize='auto', annotate=None,
            show_components=None, verbose=True, **kwargs):
        """
        Fit multiple gaussians (or other profiles)

        fittype - What function will be fit?  fittype must have been Registryed in the
            singlefitters dict.  Uses default ('gaussian') if not specified
        renormalize - if 'auto' or True, will attempt to rescale small data (<1e-9) to be 
            closer to 1 (scales by the median) so that the fit converges better
        """
        self.setfitspec()
        #if self.fitkwargs.has_key('negamp'): self.fitkwargs.pop('negamp') # We now do this in gaussfitter.py
        if fittype is not None: self.fittype = fittype
        if len(self.guesses) < self.Registry.npars[self.fittype]:
            raise ValueError("Too few parameters input.  Need at least %i for %s models" % (self.Registry.npars[self.fittype],self.fittype))
        self.npeaks = len(self.guesses)/self.Registry.npars[self.fittype]
        self.fitter = self.Registry.multifitters[self.fittype]
        self.vheight = False

        # add kwargs to fitkwargs
        self.fitkwargs.update(kwargs)
        
        scalefactor = 1.0
        if renormalize in ('auto',True):
            datarange = self.spectofit[self.xmin:self.xmax].max() - self.spectofit[self.xmin:self.xmax].min()
            if abs(datarange) < 1e-9:
                scalefactor = np.median(np.abs(self.spectofit))
                if verbose: print "Renormalizing data by factor %e to improve fitting procedure" % scalefactor
                self.spectofit /= scalefactor
                self.errspec   /= scalefactor
                for ii in xrange(self.npeaks): # assume first parameter is amplitude
                    self.guesses[self.fitter.npars*ii] /= scalefactor

        mpp,model,mpperr,chi2 = self.fitter(
                self.Spectrum.xarr[self.xmin:self.xmax], 
                self.spectofit[self.xmin:self.xmax], 
                err=self.errspec[self.xmin:self.xmax],
                npeaks=self.npeaks,
                params=self.guesses,
                **self.fitkwargs)

        self.spectofit *= scalefactor
        self.errspec   *= scalefactor

        self.mpfit_status = models.mpfit_messages[self.fitter.mp.status]

        if model is None:
            raise ValueError("Model was not set by fitter.  Examine your fitter.")
        self.chi2 = chi2
        self.dof  = self.xmax-self.xmin-self.npeaks*self.Registry.npars[self.fittype]
        self.model = model * scalefactor
        for ii in xrange(self.npeaks): # assume first parameter is amplitude
            mpp[self.fitter.npars*ii]    *= scalefactor
            mpperr[self.fitter.npars*ii] *= scalefactor
        self.modelpars = mpp.tolist()
        self.modelerrs = mpperr.tolist()
        self.parinfo = self.fitter.parinfo
        self.residuals = self.spectofit[self.xmin:self.xmax] - self.model
        if self.Spectrum.plotter.axis is not None:
            self.plot_fit(annotate=annotate, 
                    show_components=show_components, **kwargs)
                
        # Re-organize modelerrs so that any parameters that are tied to others inherit the errors of the params they are tied to
        if 'tied' in self.fitkwargs:
            for ii, element in enumerate(self.fitkwargs['tied']):
                if not element.strip(): continue
                
                i1 = element.index('[') + 1
                i2 = element.index(']')
                loc = int(element[i1:i2])
                self.modelerrs[ii] = self.modelerrs[loc]

        # make sure the full model is populated
        self._full_model()
                
    def peakbgfit(self, usemoments=True, annotate=None, vheight=True, height=0,
            negamp=None, fittype=None, renormalize='auto', 
            show_components=None, debug=False, nsigcut_moments=None, **kwargs):
        """
        Fit a single peak (plus a background)

        usemoments - The initial guess will be set by the fitter's 'moments' function
            (this overrides 'guesses')
        annotate - Make a legend?
        vheight - Fit a (constant) background as well as a peak?
        height - initial guess for background
        negamp - If True, assumes amplitude is negative.  If False, assumes positive.  If 
            None, can be either.
        fittype - What function will be fit?  fittype must have been Registryed in the
            singlefitters dict
        renormalize - if 'auto' or True, will attempt to rescale small data (<1e-9) to be 
            closer to 1 (scales by the median) so that the fit converges better
        nsigcut_moments - pass to moment guesser; can do a sigma cut for moment guessing
        """
        self.npeaks = 1
        self.auto = True
        self.setfitspec()
        if fittype is not None: self.fittype=fittype
        if usemoments: # this can be done within gaussfit but I want to save them
            # use this INDEPENDENT of fittype for now (voigt and gauss get same guesses)
            self.guesses = self.Registry.singlefitters[self.fittype].moments(
                    self.Spectrum.xarr[self.xmin:self.xmax],
                    self.spectofit[self.xmin:self.xmax], vheight=vheight,
                    negamp=negamp, nsigcut=nsigcut_moments, **kwargs)
            #if vheight is False: self.guesses = [height]+self.guesses
        else:
            if negamp: self.guesses = [height,-1,0,1]
            else:  self.guesses = [height,1,0,1]

        NP = self.Registry.singlefitters[self.fittype].default_npars
        if NP > 3:
            for ii in xrange(3,NP):
                self.guesses += [0.0]

        self.fitter = self.Registry.singlefitters[self.fittype]

        if debug: print "n(guesses): %s  Guesses: %s  vheight: %s " % (len(self.guesses),self.guesses,vheight)

        scalefactor = 1.0
        if renormalize in ('auto',True):
            datarange = self.spectofit[self.xmin:self.xmax].max() - self.spectofit[self.xmin:self.xmax].min()
            if abs(datarange) < 1e-9:
                scalefactor = np.median(np.abs(self.spectofit))
                print "Renormalizing data by factor %e to improve fitting procedure" % scalefactor
                self.spectofit /= scalefactor
                self.errspec   /= scalefactor
                self.guesses[0] /= scalefactor
                if vheight: self.guesses[1] /= scalefactor

        if debug: print "Guesses before fit: ",self.guesses
        mpp,model,mpperr,chi2 = self.fitter(
                self.Spectrum.xarr[self.xmin:self.xmax],
                self.spectofit[self.xmin:self.xmax],
                err=self.errspec[self.xmin:self.xmax],
                vheight=vheight,
                params=self.guesses,
                debug=debug,
                **self.fitkwargs)
        if debug: print "Guesses, fits after: ",self.guesses, mpp

        self.spectofit *= scalefactor
        self.errspec   *= scalefactor
        
        self.mpfit_status = models.mpfit_messages[self.fitter.mp.status]
        self.parinfo = self.fitter.parinfo

        if model is None:
            raise ValueError("Model was not set by fitter.  Examine your fitter.")
        self.chi2 = chi2
        self.dof  = self.xmax-self.xmin-self.npeaks*self.Registry.npars[self.fittype]-vheight
        self.vheight=vheight
        if vheight: 
            self.Spectrum.baseline.baselinepars = [mpp[0]*scalefactor] # first item in list form
            self.Spectrum.baseline.basespec = self.Spectrum.data*0 + mpp[0]*scalefactor
            self.model = model*scalefactor - mpp[0]*scalefactor
            # I removed this recently for some reason, but more code depends on it being in place
            # Need to figure out *WHY* anything would want an extra parameter
            if len(mpp) == self.fitter.npars+1:
                mpp = mpp[1:]
        else: self.model = model*scalefactor
        self.residuals = self.spectofit[self.xmin:self.xmax] - self.model*scalefactor
        self.modelpars = mpp.tolist()
        self.modelerrs = mpperr.tolist()
        # ONLY the amplitude was changed
        self.modelpars[0] *= scalefactor
        self.modelerrs[0] *= scalefactor
        # this was for height, but that's now dealt with above
        #self.modelpars[1] *= scalefactor
        #self.modelerrs[1] *= scalefactor
        if self.Spectrum.plotter.axis is not None:
            self.plot_fit(annotate=annotate,
                show_components=show_components, **kwargs)

        # make sure the full model is populated
        self._full_model(debug=debug)

        if debug: print "Guesses, fits after vheight removal: ",self.guesses, mpp

    def _full_model(self, debug=False, **kwargs):
        """
        Compute the model for the whole spectrum
        """
        # requires self.modelpars to be a list... hope it is
        if self.vheight:
            if self.Spectrum.baseline.baselinepars is not None:
                mpp = [self.Spectrum.baseline.baselinepars[0]] + self.modelpars
            else:
                mpp = [0] + self.modelpars
        else:
            mpp = self.modelpars
        if debug: print "_full_model mpp: ",mpp
        self.fullmodel = self.fitter.n_modelfunc(mpp,**self.fitter.modelfunc_kwargs)(self.Spectrum.xarr)
                
    def plot_fit(self, annotate=None, show_components=None, 
        composite_fit_color='red', component_fit_color='blue', lw=0.5,
        composite_lw=0.75, component_lw=0.75, component_kwargs={},
        **kwargs):
        """
        Plot the fit.  Must have fitted something before calling this!  
        
        It will be automatically called whenever a spectrum is fit (assuming an
        axis for plotting exists)

        kwargs are passed to the fitter's components attribute
        """
        if self.Spectrum.baseline.subtracted is False and self.Spectrum.baseline.basespec is not None:
            # don't display baseline if it's included in the fit
            plot_offset = self.Spectrum.plotter.offset+(self.Spectrum.baseline.basespec * (True-self.vheight))
        else:
            plot_offset = self.Spectrum.plotter.offset

        self._full_model()
        self.modelplot += self.Spectrum.plotter.axis.plot(
                self.Spectrum.xarr,
                self.fullmodel + plot_offset,
                color=composite_fit_color, linewidth=lw)
        
        # Plot components
        if show_components:
            self.modelcomponents = self.fitter.components(self.Spectrum.xarr,self.modelpars, **component_kwargs)
            for data in self.modelcomponents:
                # can have multidimensional components
                if len(data.shape) > 1:
                    for d in data:
                        self._plotted_components += self.Spectrum.plotter.axis.plot(self.Spectrum.xarr,
                            d + plot_offset,
                            color=component_fit_color, linewidth=component_lw)
                else:
                    self._plotted_components += self.Spectrum.plotter.axis.plot(self.Spectrum.xarr,
                        data + plot_offset,
                        color=component_fit_color, linewidth=component_lw)
                
        self.Spectrum.plotter.reset_limits(**self.Spectrum.plotter.plotkwargs)
        if self.Spectrum.plotter.autorefresh:
            self.Spectrum.plotter.refresh()

        if annotate or ((annotate is None) and self.autoannotate):
            self.annotate()
            if self.vheight: self.Spectrum.baseline.annotate()

    def fullsizemodel(self):
        """
        If the gaussian was fit to a sub-region of the spectrum,
        expand it (with zeros) to fill the spectrum.  
        """

        if self.model.shape != self.Spectrum.data.shape:
            temp = np.zeros(self.Spectrum.data.shape)
            temp[self.xmin:self.xmax] = self.model
            self.model = temp
            self.residuals = self.spectofit - self.model
            self.selectregion(reset=True)

    def plotresiduals(self,fig=2,axis=None,clear=True,**kwargs):
        """
        Plot residuals of the fit.  Specify a figure or
        axis; defaults to figure(2).

        kwargs are passed to matplotlib plot
        """
        if axis is None:
            if isinstance(fig,int):
                fig=matplotlib.pyplot.figure(fig)
            self.residualaxis = matplotlib.pyplot.gca()
            if clear: self.residualaxis.clear()
        else:
            self.residualaxis = axis
            if clear: self.residualaxis.clear()
        self.residualplot = self.residualaxis.plot(self.Spectrum.xarr[self.xmin:self.xmax],
                self.residuals,drawstyle='steps-mid',
                linewidth=0.5, color='k', **kwargs)
        if self.Spectrum.plotter.xmin is not None and self.Spectrum.plotter.xmax is not None:
            self.residualaxis.set_xlim(self.Spectrum.plotter.xmin,self.Spectrum.plotter.xmax)
        self.residualaxis.set_xlabel(self.Spectrum.plotter.xlabel)
        self.residualaxis.set_ylabel(self.Spectrum.plotter.ylabel)
        self.residualaxis.set_title("Residuals")
        self.residualaxis.figure.canvas.draw()

    def annotate(self,loc='upper right',labelspacing=0.25, markerscale=0.01,
            borderpad=0.1, handlelength=0.1, handletextpad=0.1, frameon=False, **kwargs):
        """
        Add a legend to the plot showing the fitted parameters

        clearlegend() will remove the legend
        
        kwargs passed to legend
        """
        self._clearlegend()
        pl = matplotlib.collections.CircleCollection([0],edgecolors=['k'])
        if hasattr(self.fitter,'annotations'):
            self._annotation_labels = self.fitter.annotations()
        else:
            raise Exception("Fitter %s has no annotations." % self.fitter)

        #xtypename = units.unit_type_dict[self.Spectrum.xarr.xtype]
        xcharconv = units.CaseInsensitiveDict({'frequency':'\\nu', 'wavelength':'\\lambda', 'velocity':'v', 'pixels':'x'})
        xchar = xcharconv[self.Spectrum.xarr.xtype]
        self._annotation_labels = [L.replace('x',xchar) if L[1]=='x' else L for L in self._annotation_labels]

        self.fitleg = self.Spectrum.plotter.axis.legend(
                tuple([pl]*len(self._annotation_labels)),
                self._annotation_labels, loc=loc, markerscale=markerscale,
                borderpad=borderpad, handlelength=handlelength,
                handletextpad=handletextpad, labelspacing=labelspacing,
                frameon=frameon, **kwargs)
        self.Spectrum.plotter.axis.add_artist(self.fitleg)
        self.fitleg.draggable(True)
        if self.Spectrum.plotter.autorefresh: self.Spectrum.plotter.refresh()

    def print_fit(self, print_baseline=True, **kwargs):
        """
        Print the best-fit parameters to the command line
        """

        if self.Spectrum.baseline.baselinepars is not None and print_baseline:
            print "Baseline: " + " + ".join(["%12g x^%i" % (x,i) for i,x in enumerate(self.Spectrum.baseline.baselinepars[::-1])])

        for i,p in enumerate(self.parinfo):
            print "%15s: %12g +/- %12g" % (p['parname'],p['value'],p['error'])

    def clear(self, legend=True, components=True):
        """
        Remove the fitted model from the plot

        Also removes the legend by default
        """
        if self.Spectrum.plotter.axis is not None:
            for p in self.modelplot:
                p.set_visible(False)
            if legend: self._clearlegend()
            if components: self._clearcomponents()
            if self.Spectrum.plotter.autorefresh: self.Spectrum.plotter.refresh()

    def _clearcomponents(self):
        for pc in self._plotted_components:
            pc.set_visible(False)
            if pc in self.Spectrum.plotter.axis.lines:
                self.Spectrum.plotter.axis.lines.remove(pc)
        if self.Spectrum.plotter.autorefresh: self.Spectrum.plotter.refresh()

    def _clearlegend(self):
        """
        Remove the legend from the plot window
        """
        if self.fitleg is not None: 
            # don't remove fitleg unless it's in the current axis self.fitleg.set_visible(False)
            if self.fitleg in self.Spectrum.plotter.axis.artists:
                self.Spectrum.plotter.axis.artists.remove(self.fitleg)
        if self.Spectrum.plotter.autorefresh: self.Spectrum.plotter.refresh()
    

    def savefit(self):
        """
        Save the fit parameters from a Gaussian fit to the FITS header
        ***THESE SHOULD BE WRITTEN FOR EACH TYPE OF MODEL TO BE FIT***
        """
        if self.modelpars is not None and hasattr(self.Spectrum,'header'):
            for ii,p in enumerate(self.modelpars):
                if ii % 3 == 0: self.Spectrum.header.update('AMP%1i' % (ii/3),p,comment="Gaussian best fit amplitude #%i" % (ii/3))
                if ii % 3 == 1: self.Spectrum.header.update('CEN%1i' % (ii/3),p,comment="Gaussian best fit center #%i" % (ii/3))
                if ii % 3 == 2: self.Spectrum.header.update('WID%1i' % (ii/3),p,comment="Gaussian best fit width #%i" % (ii/3))

    def downsample(self,factor):
        """
        Downsample the model spectrum (and the spectofit spectra)
        This should only be done when Spectrum.smooth is called
        """
        if self.model is not None:
            self.model = self.model[::factor]
            self.residuals = self.residuals[::factor]
        self.spectofit = self.spectofit[::factor]
        self.errspec = self.errspec[::factor]

    def crop(self,x1pix,x2pix):
        """
        When spectrum.crop is called, this must be too
        """
        if self.model is not None:
            self.model = self.model[x1pix:x2pix]

    def integral(self, direct=False, threshold='auto', integration_limits=[],
            return_error=False, **kwargs):
        """
        Return the integral of the fitted spectrum

        if direct=True, return the integral of the spectrum over a range
        defined by the threshold or integration limits if defined

        note that integration_limits will operate directly on the DATA, which means that
        if you've baselined without subtract=True, the baseline will be included in the integral

        if return_error is set, the error computed by
        sigma = sqrt(sum(sigma_i^2)) * dx
        will be returned as well
        """

        dx = self.Spectrum.xarr.cdelt()
        if dx is None:
            dx = np.abs(np.concatenate([np.diff(self.Spectrum.xarr),[0]]))
            warn("Irregular X-axis.  The last pixel is ignored.")
        else:
            # shouldn't shape be a 'propery'
            dx = np.repeat(np.abs(dx), self.Spectrum.shape())

        if direct:
            self._full_model()
            if len(integration_limits) == 2:
                x1 = np.argmin(np.abs(integration_limits[0]-self.Spectrum.xarr))
                x2 = np.argmin(np.abs(integration_limits[1]-self.Spectrum.xarr))
                if x1>x2: x1,x2 = x2,x1
                integ = (self.Spectrum.data[x1:x2] * dx[x1:x2]).sum()
                if return_error:
                    # compute error assuming a "known mean" (not a sample mean).  If sample mean, multiply
                    # by sqrt(len(dx)/(len(dx)-1))  (which should be very near 1)
                    error = np.sqrt((dx[x1:x2] * self.Spectrum.error[x1:x2]**2).sum() / dx[x1:x2].sum())
                    return integ,error
                else:
                    return integ
            elif threshold=='auto':
                threshold = 0.01 * np.abs( self.fullmodel ).max()

            OK = np.abs( self.fullmodel ) > threshold
            integ = (self.spectofit[OK] * dx[OK]).sum()
            error = np.sqrt((self.errspec[OK]**2 * dx[OK]).sum()/dx[OK].sum())
        else:
            if not hasattr(self.fitter,'integral'):
                raise AttributeError("The fitter %s does not have an integral implemented" % self.fittype)

            if self.Spectrum.xarr.cdelt() is not None:
                dx = np.median(dx)
                integ = self.fitter.integral(self.modelpars, **kwargs) * dx
                if return_error:
                    if mycfg.WARN: print "WARNING: The computation of the error on the integral is not obviously correct or robust... it's just a guess."
                    self._full_model()
                    OK = np.abs( self.fullmodel ) > threshold
                    error = np.sqrt((self.errspec[OK]**2).sum()) * dx
                    #raise NotImplementedError("We haven't written up correct error estimation for integrals of fits")
            else:
                integ = 0
                error = 0
                warn("An analytic integal could not be computed.  Try direct=True when integrating, or find a way to linearize the X-axis")
        if return_error:
            return integ,error
        else:
            return integ

    def shift_pars(self, frame=None):
        """
        Shift the velocity / wavelength / frequency of the fitted parameters
        into a different frame

        Right now this only takes care of redshift and only if redshift is defined.
        It should be extended to do other things later
        """
        for ii,pi in enumerate(self.parinfo):
            for partype in ('shift','offset','velo'):
                if partype in str.lower(pi['parname']):
                    if frame is not None:
                        self.modelpars[ii] = self.Spectrum.xarr.x_in_frame(self.modelpars[ii], frame)

    def moments(self, **kwargs):
        """
        Return the moments 

        see :mod:`moments`
        """
        return self.Registry.singlefitters[self.fittype].moments(
                self.Spectrum.xarr[self.xmin:self.xmax],
                self.spectofit[self.xmin:self.xmax],  **kwargs)

    def button3action(self, event, debug=False):
        """
        Disconnect the interactiveness
        Perform the fit (or die trying)
        Hide the guesses
        """
        self.Spectrum.plotter.figure.canvas.mpl_disconnect(self.click)
        self.Spectrum.plotter.figure.canvas.mpl_disconnect(self.keyclick)
        if self.npeaks > 0:
            print len(self.guesses)/3," Guesses: ",self.guesses," X channel range: ",self.xmin,self.xmax
            if len(self.guesses) % 3 == 0:
                self.multifit()
                for p in self.button2plot + self.button1plot:
                    p.set_visible(False)
            else: 
                print "error, wrong # of pars"

    def copy(self, parent=None):
        """
        Create a copy of the spectral fit - includes copies of the _full_model,
        the registry, the fitter, parinfo, modelpars, modelerrs, model, npeaks

        [ parent ] 
            A spectroscopic axis instance that is the parent of the specfit
            instance.  This needs to be specified at some point, but defaults
            to None to prevent overwriting a previous plot.
        """

        newspecfit = copy.copy(self)
        newspecfit.Spectrum = parent
        newspecfit.modelpars = self.modelpars
        newspecfit.modelerrs = self.modelerrs
        newspecfit.model = self.model
        newspecfit.npeaks = self.npeaks
        if hasattr(self,'parinfo'):
            newspecfit.parinfo = copy.copy( self.parinfo )
        if hasattr(self,'fitter'):
            newspecfit.fitter = copy.copy( self.fitter )
        if hasattr(self,'fullmodel'):
            newspecfit._full_model()

        if parent is not None:
            newspecfit.Spectrum.plotter = parent.plotter
        else:
            newspecfit.Spectrum.plotter = None

        return newspecfit
