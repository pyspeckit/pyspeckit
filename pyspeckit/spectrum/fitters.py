import matplotlib
import matplotlib.cbook as mpcb
import matplotlib.pyplot as pyplot
import numpy as np
from config import spcfg

class Registry(object):
    """
    This class is a simple wrapper to prevent fitter properties from being globals
    """

    def __init__(self):
        self.npars = {}
        self.multifitters = {}
        self.singlefitters = {}
        self.writers = {}
        self.fitkeys = {}

        self.interactive_help_message = """
        Left-click or hit 'p' twice to select a fitting range, then middle-click or hit
        'm' twice to select a peak and width.  When you're done, right-click or hit 'd'
        to perform the fit and disconnect the mouse and keyboard.  '?' will print this
        help message again.
        You can select different fitters to use with the interactive fitting routine.
        The default is gaussian ('g')
        """


class Specfit(object):

    def __init__(self, Spectrum, autoannotate=bool(spcfg.cfg['annotate']),
            Registry=None):
        self.model = None
        self.modelpars = None
        self.modelerrs = None
        self.modelplot = None
        self.modelcomponents = None
        self.guessplot = []
        self.fitregion = []
        self.npeaks = 0
        self.nclicks_b1 = 0
        self.nclicks_b2 = 0
        self.gx1 = 0
        self.gx2 = Spectrum.data.shape[0]
        self.guesses = []
        self.click = 0
        self.fitkwargs = {}
        self.auto = False
        self.Spectrum = Spectrum
        self.specplotter = self.Spectrum.plotter
        self.fitleg=None
        self.residuals=None
        self.setfitspec()
        self.fittype = 'gaussian'
        self.measurements = None
        self.vheight=None
        self.Registry = Registry
        #self.seterrspec()
        
        # config file stuff
        self.cfg = spcfg.cfg
        self.fitcolor = self.cfg['fit_color']
        self.compcolor = self.cfg['comp_color']
        self.fitlw = self.cfg['fit_lw']
        self.complw = self.cfg['comp_lw']
        self.autoannotate = autoannotate
        self.show_components = bool(self.cfg['show_components'])

    def __call__(self, interactive=False, usemoments=True, fitcolor=None,
            multifit=False, guesses=None, annotate=None, save=True,
            fittype='gaussian', compcolor = None, fitlw = None, complw = None, 
            debug=False, clear_all_connections=False, **kwargs):
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
  
        # don't change defaults if fitcolor is not None: self.fitcolor = fitcolor
        # don't change defaults if compcolor is not None: self.compcolor = compcolor
        # don't change defaults if fitlw is not None: self.fitlw = fitlw
        # don't change defaults if complw is not None: self.complw = complw
        if annotate is not None: self.autoannotate = annotate

        self.clear()
        self.selectregion(**kwargs)
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
            if self.specplotter.axis is None:
                raise Exception("Interactive fitting requires a plotter.")
            else:
                self.specplotter.axis.set_autoscale_on(False)
            print self.Registry.interactive_help_message
            self.nclicks_b1 = 0
            self.nclicks_b2 = 0
            self.guesses = []
            makeguess = self.makeguess_debug if debug else self.makeguess
            
            if clear_all_connections:
                # this is really ugly, but needs to be done in order to prevent multiple overlapping calls...
                cids_to_remove = []
                for eventtype in ('button_press_event','key_press_event'):
                    for key,val in self.specplotter.figure.canvas.callbacks.callbacks[eventtype].iteritems():
                        if "makeguess" in val.func.__name__:
                            cids_to_remove.append(key)
                            if debug: print "Removing CID #%i with attached function %s" % (key,val.func.__name__)
                for cid in cids_to_remove:
                    self.specplotter.figure.canvas.mpl_disconnect(cid)

            self.click = self.specplotter.axis.figure.canvas.mpl_connect('button_press_event',makeguess)
            self.keyclick = self.specplotter.axis.figure.canvas.mpl_connect('key_press_event',makeguess)
        elif multifit and self.fittype in self.Registry.multifitters or guesses is not None:
            if guesses is None:
                print "You must input guesses when using multifit.  Also, baseline (continuum fit) first!"
                return
            else:
                self.guesses = guesses
                self.multifit()
        # SINGLEFITTERS SHOULD BE PHASED OUT
        elif self.fittype in self.Registry.singlefitters:
            #print "Non-interactive, 1D fit with automatic guessing"
            if self.Spectrum.baseline.order is None:
                self.Spectrum.baseline.order=0
                self.peakbgfit(usemoments=usemoments, **kwargs)
            else:
                self.peakbgfit(usemoments=usemoments, 
                        vheight=False, height=0.0, **kwargs)
            if self.specplotter.autorefresh: self.specplotter.refresh()
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
            print "WARNING: Baseline / continuum is negative: equivalent width is poorly defined."
        diffspec = (self.Spectrum.baseline.basespec - self.Spectrum.data)
        dx = np.abs((self.Spectrum.xarr[self.gx2-1]-self.Spectrum.xarr[self.gx1]) / (self.gx2-self.gx1))
        sumofspec = diffspec[self.gx1:self.gx2].sum() * dx
        eqw = sumofspec / np.median(self.Spectrum.baseline.basespec)
        if plot:
            midpt_pixel = np.round((self.gx1+self.gx2)/2.0)
            midpt       = self.Spectrum.xarr[midpt_pixel]
            midpt_level = self.Spectrum.baseline.basespec[midpt_pixel]
            print "EQW plotting: ",midpt,midpt_pixel,midpt_level,eqw
            self.specplotter.axis.fill_between(
                    [midpt-eqw/2.0,midpt+eqw/2.0],
                    [0,0],
                    [midpt_level,midpt_level],
                    color=plotcolor,
                    alpha=alpha,
                    label='EQW: %0.3g' % eqw)
            if annotate:
                self.specplotter.axis.legend(
                        [(matplotlib.collections.CircleCollection([0],facecolors=[plotcolor],edgecolors=[plotcolor]))],
                        [('EQW: %0.3g' % eqw)], 
                        markerscale=0.01, borderpad=0.1, handlelength=0.1,
                        handletextpad=0.1, loc=loc)
            if self.specplotter.autorefresh:
                self.specplotter.refresh()
        return eqw

    
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
            if self.Spectrum.baseline.subtracted is False and self.Spectrum.baseline.basespec is not None:
                self.spectofit -= self.Spectrum.baseline.basespec
        OKmask = (self.spectofit==self.spectofit)
        self.spectofit[(True-OKmask)] = 0
        self.seterrspec()
        self.errspec[(True-OKmask)] = 1e10

    def multifit(self, fittype=None, renormalize='auto', annotate=None):
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
        self.npeaks = len(self.guesses)/self.Registry.npars[self.fittype]
        self.fitter = self.Registry.multifitters[self.fittype]
        self.vheight = False
        
        scalefactor = 1.0
        if renormalize in ('auto',True):
            datarange = self.spectofit[self.gx1:self.gx2].max() - self.spectofit[self.gx1:self.gx2].min()
            if abs(datarange) < 1e-9:
                scalefactor = np.median(np.abs(self.spectofit))
                print "Renormalizing data by factor %e to improve fitting procedure" % scalefactor
                self.spectofit /= scalefactor
                self.errspec   /= scalefactor
                for ii in xrange(self.npeaks): # assume first parameter is amplitude
                    self.guesses[self.fitter.npars*ii] /= scalefactor

        mpp,model,mpperr,chi2 = self.fitter(
                self.Spectrum.xarr[self.gx1:self.gx2], 
                self.spectofit[self.gx1:self.gx2], 
                err=self.errspec[self.gx1:self.gx2],
                npeaks=self.npeaks,
                params=self.guesses,
                **self.fitkwargs)

        self.spectofit *= scalefactor
        self.errspec   *= scalefactor

        if model is None:
            raise ValueError("Model was not set by fitter.  Examine your fitter.")
        self.chi2 = chi2
        self.dof  = self.gx2-self.gx1-self.npeaks*self.Registry.npars[self.fittype]
        self.model = model * scalefactor
        for ii in xrange(self.npeaks): # assume first parameter is amplitude
            mpp[self.fitter.npars*ii]    *= scalefactor
            mpperr[self.fitter.npars*ii] *= scalefactor
        self.modelpars = mpp.tolist()
        self.modelerrs = mpperr.tolist()
        self.residuals = self.spectofit[self.gx1:self.gx2] - self.model
        if self.specplotter.axis is not None:
            self.plot_fit(annotate=annotate)
                
        # Re-organize modelerrs so that any parameters that are tied to others inherit the errors of the params they are tied to
        if self.fitkwargs.has_key('tied'):
            for ii, element in enumerate(self.fitkwargs['tied']):
                if not element.strip(): continue
                
                i1 = element.index('[') + 1
                i2 = element.index(']')
                loc = int(element[i1:i2])
                self.modelerrs[ii] = self.modelerrs[loc]
                
    
    def peakbgfit(self, usemoments=True, annotate=None, vheight=True, height=0,
            negamp=None, fittype=None, renormalize='auto', **kwargs):
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
        """
        self.npeaks = 1
        self.auto = True
        self.setfitspec()
        if fittype is not None: self.fittype=fittype
        if usemoments: # this can be done within gaussfit but I want to save them
            # use this INDEPENDENT of fittype for now (voigt and gauss get same guesses)
            self.guesses = self.Registry.singlefitters[self.fittype].moments(
                    self.Spectrum.xarr[self.gx1:self.gx2],
                    self.spectofit[self.gx1:self.gx2],
                    vheight=vheight,negamp=negamp,**kwargs)
            if vheight is False: self.guesses = [height]+self.guesses
        else:
            if negamp: self.guesses = [height,-1,0,1]
            else:  self.guesses = [height,1,0,1]
        if self.fittype == 'voigt':
            self.guesses += [0.0]
        self.fitter = self.Registry.singlefitters[self.fittype]

        scalefactor = 1.0
        if renormalize in ('auto',True):
            datarange = self.spectofit[self.gx1:self.gx2].max() - self.spectofit[self.gx1:self.gx2].min()
            if abs(datarange) < 1e-9:
                scalefactor = np.median(np.abs(self.spectofit))
                print "Renormalizing data by factor %e to improve fitting procedure" % scalefactor
                self.spectofit /= scalefactor
                self.errspec   /= scalefactor
                self.guesses[0] /= scalefactor
                self.guesses[1] /= scalefactor

        mpp,model,mpperr,chi2 = self.fitter(
                self.Spectrum.xarr[self.gx1:self.gx2],
                self.spectofit[self.gx1:self.gx2],
                err=self.errspec[self.gx1:self.gx2],
                vheight=vheight,
                params=self.guesses,
                **self.fitkwargs)

        self.spectofit *= scalefactor
        self.errspec   *= scalefactor

        if model is None:
            raise ValueError("Model was not set by fitter.  Examine your fitter.")
        self.chi2 = chi2
        self.dof  = self.gx2-self.gx1-self.npeaks*self.Registry.npars[self.fittype]-vheight
        if vheight: 
            self.vheight=True
            self.Spectrum.baseline.baselinepars = mpp[:1]*scalefactor # first item in list form
            self.Spectrum.baseline.basespec = self.Spectrum.data*0 + mpp[0]*scalefactor
            self.model = model*scalefactor - mpp[0]*scalefactor
        else: self.model = model*scalefactor
        self.residuals = self.spectofit[self.gx1:self.gx2] - self.model*scalefactor
        self.modelpars = mpp[1:].tolist()
        self.modelerrs = mpperr[1:].tolist()
        self.modelpars[0] *= scalefactor
        self.modelerrs[0] *= scalefactor
        self.modelpars[1] *= scalefactor
        self.modelerrs[1] *= scalefactor
        if self.specplotter.axis is not None:
            self.plot_fit(annotate=annotate)

    def plot_fit(self, annotate=None, show_components=None):
        """
        Plot the fit.  Must have fitted something before calling this!  
        
        It will be automatically called whenever a spectrum is fit (assuming an
        axis for plotting exists)
        """
        if self.Spectrum.baseline.subtracted is False and self.Spectrum.baseline.basespec is not None:
            plotmodel = self.model+self.specplotter.offset+self.Spectrum.baseline.basespec[self.gx1:self.gx2]
        else:
            plotmodel = self.model+self.specplotter.offset
        self.modelplot = self.specplotter.axis.plot(
                self.Spectrum.xarr[self.gx1:self.gx2],
                plotmodel,
                color=self.fitcolor, linewidth=self.fitlw)
        
        if show_components is not None:
            self.show_components = show_components
        # Plot components
        if self.show_components:
            self.modelcomponents = self.fitter.components(self.Spectrum.xarr[self.gx1:self.gx2],self.modelpars)
            for data in self.modelcomponents:
                self.specplotter.axis.plot(self.Spectrum.xarr[self.gx1:self.gx2],
                    data+self.specplotter.offset+self.Spectrum.baseline.basespec[self.gx1:self.gx2],
                    color=self.compcolor, linewidth=self.complw)                
                
        self.specplotter.reset_limits(**self.specplotter.plotkwargs)
        if self.specplotter.autorefresh:
            self.specplotter.refresh()

        if annotate is not None:
            self.autoannotate = annotate

        if self.autoannotate:
            self.annotate()
            if self.vheight: self.Spectrum.baseline.annotate()

    def fullsizemodel(self):
        """
        If the gaussian was fit to a sub-region of the spectrum,
        expand it (with zeros) to fill the spectrum.  You can 
        always recover the original by:
        origmodel = model[gx1:gx2]
        """

        if self.model.shape != self.Spectrum.data.shape:
            temp = np.zeros(self.Spectrum.data.shape)
            temp[self.gx1:self.gx2] = self.model
            self.model = temp
            self.residuals = self.spectofit - self.model

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
        self.residualplot = self.residualaxis.plot(self.Spectrum.xarr[self.gx1:self.gx2],
                self.residuals,drawstyle='steps-mid',
                linewidth=0.5, color='k', **kwargs)
        if self.specplotter.xmin is not None and self.specplotter.xmax is not None:
            self.residualaxis.set_xlim(self.specplotter.xmin,self.specplotter.xmax)
        self.residualaxis.set_xlabel(self.specplotter.xlabel)
        self.residualaxis.set_ylabel(self.specplotter.ylabel)
        self.residualaxis.set_title("Residuals")
        self.residualaxis.figure.canvas.draw()

    def annotate(self,loc='upper right'):
        """
        Add a legend to the plot showing the fitted parameters

        clearlegend() will remove the legend
        """
        self.clearlegend()
        pl = matplotlib.collections.CircleCollection([0],edgecolors=['k'])
        if hasattr(self.fitter,'annotations'):
            labels = self.fitter.annotations()
        else:
            raise Exception("Fitter %s has no annotations." % self.fitter)
        self.fitleg = self.specplotter.axis.legend(
                tuple([pl]*self.Registry.npars[self.fittype]*self.npeaks),
                labels,
                loc=loc,markerscale=0.01,
                borderpad=0.1, handlelength=0.1, handletextpad=0.1
                )
        self.specplotter.axis.add_artist(self.fitleg)
        self.fitleg.draggable(True)
        if self.specplotter.autorefresh: self.specplotter.refresh()

    def selectregion(self,xmin=None,xmax=None,xtype='wcs', reset=False,
            debug=False, **kwargs):
        """
        Pick a fitting region in either WCS units or pixel units

        reset - if true, overrides input xmin,xmax and selects the full range
        """
        if xmin is not None and xmax is not None:
            if xtype in ('wcs','WCS','velo','velocity','wavelength','frequency','freq','wav'):
                self.gx1 = np.argmin(abs(xmin-self.Spectrum.xarr))
                self.gx2 = np.argmin(abs(xmax-self.Spectrum.xarr))
            else:
                self.gx1 = xmin
                self.gx2 = xmax
        elif self.specplotter.xmin is not None and self.specplotter.xmax is not None:
            self.gx1 = np.argmin(abs(self.specplotter.xmin-self.Spectrum.xarr))
            self.gx2 = np.argmin(abs(self.specplotter.xmax-self.Spectrum.xarr))
        elif reset:
            self.gx1 = 0
            self.gx2 = self.Spectrum.data.shape[0]
            #raise ValueError("Need to input xmin and xmax, or have them set by plotter, for selectregion.")
        else:
            print "Left region selection unchanged.  xminpix, xmaxpix: %i,%i" % (self.gx1,self.gx2)

        if self.gx1 == self.gx2:
            # Reset if there is no fitting region
            self.gx1 = 0
            self.gx2 = self.Spectrum.data.shape[0]
            if debug: print "Reset to full range because the endpoints were equal"
        elif self.gx1>self.gx2: 
            # Swap endpoints if the axis has a negative delta-X
            self.gx1,self.gx2 = self.gx2,self.gx1
            if debug: print "Swapped endpoints because the left end was greater than the right"


    def selectregion_interactive(self,event,debug=False):
        """
        ***For window-interactive use only!***
        (i.e., you probably shouldn't call this from the command line or a script)

        Defines the fitting region of the spectrum
        """
        if self.nclicks_b1 == 0:
            self.gx1 = np.argmin(abs(event.xdata-self.Spectrum.xarr))
            if debug: print "Left at %g.  Click #%i" % (self.gx1,self.nclicks_b1)
            self.nclicks_b1 += 1
        elif self.nclicks_b1 == 1:
            self.gx2 = np.argmin(abs(event.xdata-self.Spectrum.xarr))
            if debug: print "Right at %g.  Click #%i" % (self.gx2,self.nclicks_b2)
            self.nclicks_b1 -= 1
            if self.gx1 > self.gx2: self.gx1,self.gx2 = self.gx2,self.gx1
            if abs(self.gx1-self.gx2) > 3: # can't fit w/ fewer data than pars
                self.fitregion = self.specplotter.axis.plot(
                        self.Spectrum.xarr[self.gx1:self.gx2],
                        self.Spectrum.data[self.gx1:self.gx2]+self.specplotter.offset,
                        drawstyle='steps-mid',
                        color='c')
                #self.specplotter.plot(**self.specplotter.plotkwargs)
                if self.guesses == []:
                    self.guesses = self.Registry.singlefitters['gaussian'].moments(
                            self.Spectrum.xarr[self.gx1:self.gx2],
                            self.spectofit[self.gx1:self.gx2],
                            vheight=0)
                    self.npeaks = 1
                    self.auto = True
            else:
                print "Fitting region is too small (channels %i:%i).  Try again." % (self.gx1,self.gx2)

    def guesspeakwidth(self,event,debug=False):
        """
        Interactively guess the peak height and width from user input

        Width is assumed to be half-width-half-max
        """
        if self.nclicks_b2 % 2 == 0:
            if self.auto:
                self.guesses[:2] = [event.ydata,event.xdata]
            else:
                self.guesses += [event.ydata,event.xdata,1]
                self.npeaks += 1
            self.nclicks_b2 += 1
            if debug: print "Peak %i click %i at x,y %g,%g" % (self.npeaks,self.nclicks_b2,event.xdata,event.ydata)
            self.guessplot += [self.specplotter.axis.scatter(event.xdata,event.ydata,marker='x',c='r')]
            #self.specplotter.refresh() #plot(**self.specplotter.plotkwargs)
        elif self.nclicks_b2 % 2 == 1:
            self.guesses[-1] = abs(event.xdata-self.guesses[-2]) / np.sqrt(2*np.log(2))
            self.nclicks_b2 += 1
            if debug: print "Width %i click %i at x,y %g,%g" % (self.npeaks,self.nclicks_b2,event.xdata,event.ydata)
            self.guessplot += self.specplotter.axis.plot([event.xdata,
                2*self.guesses[-2]-event.xdata],[event.ydata]*2,
                color='r')
            #self.specplotter.refresh() #plot(**self.specplotter.plotkwargs)
            if self.auto:
                self.auto = False
            if self.nclicks_b2 / 2 > self.npeaks:
                print "There have been %i middle-clicks but there are only %i gaussians" % (self.nclicks_b2,self.npeaks)
                self.npeaks += 1

    def clear(self,legend=True):
        """
        Remove the fitted model from the plot

        Also removes the legend by default
        """
        if self.modelplot is not None:
            for p in self.modelplot:
                p.set_visible(False)
        if legend: self.clearlegend()

    def makeguess_debug(self,event):
        return self.makeguess(event,debug=True)

    def makeguess(self,event,debug=False):
        """
        ***For window-interactive use only!***
        (i.e., you probably shouldn't call this from the command line or a script)

        Given a set of clicks or button presses, sets the fit guesses
        """
        toolbar = self.specplotter.figure.canvas.manager.toolbar
        if toolbar.mode == '' and self.specplotter.axis in event.canvas.figure.axes:
            if hasattr(event,'button'):
                button = event.button
            elif hasattr(event,'key'):
                button = event.key

            if debug:
                print "button: ",button

            if button in ('p','P','1',1):
                self.selectregion_interactive(event,debug=debug)
            elif button in ('m','M','2',2):
                self.guesspeakwidth(event,debug=debug)
            elif button in ('d','D','3',3):
                self.specplotter.figure.canvas.mpl_disconnect(self.click)
                self.specplotter.figure.canvas.mpl_disconnect(self.keyclick)
                if self.npeaks > 0:
                    print len(self.guesses)/3," Guesses: ",self.guesses," X channel range: ",self.gx1,self.gx2
                    if len(self.guesses) % 3 == 0:
                        self.multifit()
                        for p in self.guessplot + self.fitregion:
                            p.set_visible(False)
                    else: 
                        print "error, wrong # of pars"
            elif button in ('?'):
                print self.Registry.interactive_help_message
            elif button in self.Registry.fitkeys:
                fittername = self.Registry.fitkeys[button]
                if fittername in self.Registry.multifitters:
                    self.fitter = self.Registry.multifitters[fittername]
                    self.fittype = fittername
                    print "Selected multi-fitter %s" % fittername
                elif fittername in self.Registry.singlefitters:
                    self.fitter = self.Registry.singlefitters[fittername]
                    self.fittype = fittername
                    print "Selected single-fitter %s" % fittername
                else: 
                    print "ERROR: Did not find fitter %s" % fittername
            if self.specplotter.autorefresh: self.specplotter.refresh()

    def clearlegend(self):
        """
        Remove the legend from the plot window
        """
        if self.fitleg is not None: 
            self.fitleg.set_visible(False)
            if self.fitleg in self.specplotter.axis.artists:
                self.specplotter.axis.artists.remove(self.fitleg)
        if self.specplotter.autorefresh: self.specplotter.refresh()
    
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

        if not hasattr(self.fitter,'integral'):
            raise AttributeError("The fitter %s does not have an integral implemented" % self.fittype)

        if direct:
            dx = np.abs(self.Spectrum.xarr.cdelt())
            if len(integration_limits) == 2:
                x1 = np.argmin(np.abs(integration_limits[0]-self.Spectrum.xarr))
                x2 = np.argmin(np.abs(integration_limits[1]-self.Spectrum.xarr))
                if x1>x2: x1,x2 = x2,x1
                integ = self.Spectrum.data[x1:x2] * dx
                if return_error:
                    error = np.sqrt((self.Spectrum.error[x1:x2]**2).sum()) * dx
                    return integ,error
                else:
                    return integ
            elif threshold=='auto':
                threshold = 0.01 * np.abs( self.model ).max()

            OK = np.abs( self.model ) > threshold
            integ = self.spectofit[OK].sum() * dx
            error = np.sqrt((self.errspec[OK]**2).sum()) * dx
        else:
            integ = self.fitter.integral(self.modelpars, **kwargs)
            if return_error:
                raise NotImplementedError("We haven't written up correct error estimation for integrals of fits")
        if return_error:
            return integ,error
        else:
            return integ


