import matplotlib
import matplotlib.cbook as mpcb
import matplotlib.pyplot as pyplot
import numpy as np
from config import *

interactive_help_message = """
Left-click or hit 'p' twice to select a fitting range, then middle-click or hit
'm' twice to select a peak and width.  When you're done, right-click or hit 'd'
to perform the fit and disconnect the mouse and keyboard.  '?' will print this
help message again.
"""

npars = {}
multifitters = {}
singlefitters = {}
writers = {}

class Specfit(object):

    def __init__(self,Spectrum):
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
        self.fitxarr = None
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
        #self.seterrspec()
        
        # config file stuff
        self.cfg = spcfg.cfg
        self.fitcolor = self.cfg['fit_color']
        self.compcolor = self.cfg['comp_color']
        self.fitlw = self.cfg['fit_lw']
        self.complw = self.cfg['comp_lw']
        self.autoannotate = bool(self.cfg['annotate'])
        self.show_components = bool(self.cfg['show_components'])

    def __call__(self, interactive=False, usemoments=True, fitcolor=None,
            multifit=False, guesses=None, annotate=None, save=True,
            fittype='gaussian', compcolor = None, fitlw = None, complw = None, 
            **kwargs):
        """
        Fit gaussians to a spectrum

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
  
        if fitcolor is not None: self.fitcolor = fitcolor
        if compcolor is not None: self.compcolor = compcolor
        if fitlw is not None: self.fitlw = fitlw
        if complw is not None: self.complw = complw
        if annotate is None: annotate = self.autoannotate

        self.clear()
        self.selectregion(**kwargs)
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
            print interactive_help_message
            self.nclicks_b1 = 0
            self.nclicks_b2 = 0
            self.guesses = []
            self.click = self.specplotter.axis.figure.canvas.mpl_connect('button_press_event',self.makeguess)
            self.keyclick = self.specplotter.axis.figure.canvas.mpl_connect('key_press_event',self.makeguess)
        elif multifit and fittype in multifitters:
            if guesses is None:
                print "You must input guesses when using multifit.  Also, baseline (continuum fit) first!"
                return
            else:
                self.guesses = guesses
                self.multifit(fittype=fittype)
        elif fittype in singlefitters:
            #print "Non-interactive, 1D fit with automatic guessing"
            if self.Spectrum.baseline.order is None:
                self.Spectrum.baseline.order=0
                self.peakbgfit(usemoments=usemoments, annotate=annotate,
                        fittype=fittype, **kwargs)
            else:
                self.peakbgfit(usemoments=usemoments, annotate=annotate,
                        vheight=False, height=0.0, fittype=fittype, **kwargs)
            if self.specplotter.autorefresh: self.specplotter.refresh()
        else:
            if multifit:
                print "Can't fit with given fittype %s: it is not registered as a multifitter." % fittype
            else:
                print "Can't fit with given fittype %s: it is not registered as a singlefitter." % fittype
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

    def multifit(self, fittype='gaussian'):
        """
        Fit multiple gaussians (or other profiles)
        """
        self.npeaks = len(self.guesses)/npars[fittype]
        self.setfitspec()
        if self.fitkwargs.has_key('negamp'): self.fitkwargs.pop('negamp')
        self.fitter = multifitters[fittype]
        mpp,model,mpperr,chi2 = self.fitter(
                self.Spectrum.xarr[self.gx1:self.gx2], 
                self.spectofit[self.gx1:self.gx2], 
                err=self.errspec[self.gx1:self.gx2],
                npeaks=self.npeaks,
                params=self.guesses,
                **self.fitkwargs)
        if model is None:
            raise ValueError("Model was not set by fitter.  Examine your fitter.")
        self.chi2 = chi2
        self.dof  = self.gx2-self.gx1-self.npeaks*npars[fittype]
        self.model = model
        self.modelpars = mpp.tolist()
        self.modelerrs = mpperr.tolist()
        self.residuals = self.spectofit[self.gx1:self.gx2] - self.model
        if self.specplotter.axis is not None:
            self.plot_fit()
            if self.autoannotate:
                self.annotate()
    
    def peakbgfit(self, usemoments=True, annotate=True, vheight=True, height=0,
            negamp=None, fittype='gaussian', **kwargs):
        """
        Fit a single peak
        """
        self.npeaks = 1
        self.auto = True
        self.setfitspec()
        if usemoments: # this can be done within gaussfit but I want to save them
            # use this INDEPENDENT of fittype for now (voigt and gauss get same guesses)
            self.guesses = singlefitters[fittype].moments(
                    self.Spectrum.xarr[self.gx1:self.gx2],
                    self.spectofit[self.gx1:self.gx2],
                    vheight=vheight,negamp=negamp,**kwargs)
            if vheight is False: self.guesses = [height]+self.guesses
        else:
            if negamp: self.guesses = [height,-1,0,1]
            else:  self.guesses = [height,1,0,1]
        if fittype == 'voigt':
            self.guesses += [0.0]
        self.fitter = singlefitters[fittype]
        mpp,model,mpperr,chi2 = self.fitter(
                self.Spectrum.xarr[self.gx1:self.gx2],
                self.spectofit[self.gx1:self.gx2],
                err=self.errspec[self.gx1:self.gx2],
                vheight=vheight,
                params=self.guesses,
                **self.fitkwargs)
        if model is None:
            raise ValueError("Model was not set by fitter.  Examine your fitter.")
        self.chi2 = chi2
        self.dof  = self.gx2-self.gx1-self.npeaks*npars[fittype]-vheight
        if vheight: 
            self.Spectrum.baseline.baselinepars = mpp[:1] # first item in list form
            self.Spectrum.baseline.basespec = self.Spectrum.data*0 + mpp[0]
            self.model = model - mpp[0]
        else: self.model = model
        self.residuals = self.spectofit[self.gx1:self.gx2] - self.model
        self.modelpars = mpp[1:].tolist()
        self.modelerrs = mpperr[1:].tolist()
        if self.specplotter.axis is not None:
            self.plot_fit()
            if annotate:
                self.annotate()
                if vheight: self.Spectrum.baseline.annotate()

    def plot_fit(self):
        if self.Spectrum.baseline.subtracted is False and self.Spectrum.baseline.basespec is not None:
            plotmodel = self.model+self.specplotter.offset+self.Spectrum.baseline.basespec[self.gx1:self.gx2]
        else:
            plotmodel = self.model+self.specplotter.offset
        self.modelplot = self.specplotter.axis.plot(
                self.Spectrum.xarr[self.gx1:self.gx2],
                plotmodel,
                color=self.fitcolor, linewidth=self.fitlw)
        
        # Compute individual Gaussian components
        self.modelcomponents = np.empty(len(self.modelpars) / 3, np.ndarray)
        for i in np.arange(len(self.modelpars) / 3):
            self.modelcomponents[i] = singlefitters['gaussian'].onepeakgaussian(self.Spectrum.xarr[self.gx1:self.gx2],
            0.0,self.modelpars[3*i],self.modelpars[3*i+1],self.modelpars[3*i+2])
        
        # Plot components
        if self.show_components:
            for i in np.arange(len(self.modelpars) / 3):
                self.specplotter.axis.plot(self.Spectrum.xarr[self.gx1:self.gx2],
                self.modelcomponents[i], color=self.compcolor, linewidth=self.complw)                
                
        self.specplotter.plot(**self.specplotter.plotkwargs)

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
        self.residualaxis.figure.canvas.draw()

    def annotate(self,loc='upper right'):
        #text(xloc,yloc     ,"c=%g" % self.modelpars[1],transform = self.specplotter.axis.transAxes)
        #text(xloc,yloc-0.05,"w=%g" % self.modelpars[2],transform = self.specplotter.axis.transAxes)
        #text(xloc,yloc-0.10,"a=%g" % self.modelpars[0],transform = self.specplotter.axis.transAxes)
        self.clearlegend()
        pl = matplotlib.collections.CircleCollection([0],edgecolors=['k'])
        if hasattr(self.fitter,'annotations'):
            labels = self.fitter.annotations()
        else:
            raise Exception("Fitter %s has no annotations." % self.fitter)
        self.fitleg = self.specplotter.axis.legend(
                tuple([pl]*npars[self.fittype]*self.npeaks),
                labels,
                loc=loc,markerscale=0.01,
                borderpad=0.1, handlelength=0.1, handletextpad=0.1
                )
        self.fitleg.draggable(True)
        self.specplotter.axis.add_artist(self.fitleg)
        if self.specplotter.autorefresh: self.specplotter.refresh()

    def selectregion(self,xmin=None,xmax=None,xtype='wcs',**kwargs):
        """
        Pick a fitting region in either WCS units or pixel units
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
        else:
            self.gx1 = 0
            self.gx2 = self.Spectrum.data.shape[0]
            #raise ValueError("Need to input xmin and xmax, or have them set by plotter, for selectregion.")

        if self.gx1 == self.gx2:
            # Reset if there is no fitting region
            self.gx1 = 0
            self.gx2 = self.Spectrum.data.shape[0]
        elif self.gx1>self.gx2: 
            # Swap endpoints if the axis has a negative delta-X
            self.gx1,self.gx2 = self.gx2,self.gx1


    def selectregion_interactive(self,event):
        if self.nclicks_b1 == 0:
            self.gx1 = np.argmin(abs(event.xdata-self.Spectrum.xarr))
            self.nclicks_b1 += 1
        elif self.nclicks_b1 == 1:
            self.gx2 = np.argmin(abs(event.xdata-self.Spectrum.xarr))
            self.nclicks_b1 -= 1
            if self.gx1 > self.gx2: self.gx1,self.gx2 = self.gx2,self.gx1
            if abs(self.gx1-self.gx2) > 3: # can't fit w/ fewer data than pars
                self.fitregion = self.specplotter.axis.plot(
                        self.Spectrum.xarr[self.gx1:self.gx2],
                        self.Spectrum.data[self.gx1:self.gx2]+self.specplotter.offset,
                        drawstyle='steps-mid',
                        color='c')
                self.specplotter.plot(**self.specplotter.plotkwargs)
                if self.guesses == []:
                    self.guesses = singlefitters['gaussian'].moments(
                            self.Spectrum.xarr[self.gx1:self.gx2],
                            self.spectofit[self.gx1:self.gx2],
                            vheight=0)
                    self.npeaks = 1
                    self.auto = True
            else:
                print "Fitting region is too small (channels %i:%i).  Try again." % (self.gx1,self.gx2)

    def guesspeakwidth(self,event):
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
            self.guessplot += [self.specplotter.axis.scatter(event.xdata,event.ydata,marker='x',c='r')]
            self.specplotter.plot(**self.specplotter.plotkwargs)
        elif self.nclicks_b2 % 2 == 1:
            self.guesses[-1] = abs(event.xdata-self.guesses[-2]) / np.sqrt(2*np.log(2))
            self.nclicks_b2 += 1
            self.guessplot += self.specplotter.axis.plot([event.xdata,
                2*self.guesses[-2]-event.xdata],[event.ydata]*2,
                color='r')
            self.specplotter.plot(**self.specplotter.plotkwargs)
            if self.auto:
                self.auto = False
            if self.nclicks_b2 / 2 > self.npeaks:
                print "There have been %i middle-clicks but there are only %i gaussians" % (self.nclicks_b2,self.npeaks)
                self.npeaks += 1

    def clear(self,legend=True):
        if self.modelplot is not None:
            for p in self.modelplot:
                p.set_visible(False)
        if legend: self.clearlegend()

    def makeguess(self,event):
        if hasattr(event,'button'):
            button = event.button
        elif hasattr(event,'key'):
            button = event.key

        if button in ('p','P','1',1):
            self.selectregion_interactive(event)
        elif button in ('m','M','2',2):
            self.guesspeakwidth(event)
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
            print interactive_help_message
        if self.specplotter.autorefresh: self.specplotter.refresh()

    def clearlegend(self):
        if self.fitleg is not None: 
            self.fitleg.set_visible(False)
            if self.fitleg in self.specplotter.axis.artists:
                self.specplotter.axis.artists.remove(self.fitleg)
        if self.specplotter.autorefresh: self.specplotter.refresh()
    
    def savefit(self):
        if self.modelpars is not None and hasattr(self.Spectrum,'header'):
            for ii,p in enumerate(self.modelpars):
                if ii % 3 == 0: self.Spectrum.header.update('AMP%1i' % (ii/3),p,comment="Gaussian best fit amplitude #%i" % (ii/3))
                if ii % 3 == 1: self.Spectrum.header.update('CEN%1i' % (ii/3),p,comment="Gaussian best fit center #%i" % (ii/3))
                if ii % 3 == 2: self.Spectrum.header.update('WID%1i' % (ii/3),p,comment="Gaussian best fit width #%i" % (ii/3))
