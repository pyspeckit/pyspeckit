import numpy as np
import matplotlib

class Baseline:
    def __init__(self,Spectrum):
        self.baselinepars  = None
        self.order = None
        self.basespec = np.zeros(Spectrum.data.shape[0])
        self.excludemask = np.zeros(Spectrum.data.shape[0],dtype='bool')
        self.OKmask = np.ones(Spectrum.data.shape[0],dtype='bool')
        self.Spectrum = Spectrum
        self.specplotter = Spectrum.plotter
        self.blleg = None
        self.click = 0
        self.nclicks_b1 = 0
        self.nclicks_b2 = 0
        self.fitregion=[]
        self.excludevelo = []
        self.excludepix  = []
        self.subtracted = False

    def __call__(self, order=1, annotate=False, excludefit=False, save=True,
            exclude=None, exclusionlevel=0.01, interactive=False, **kwargs):
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

        "basespec" is added back to the spectrum before fitting so you can run this
        procedure multiple times without losing information
        """
        specfit = self.Spectrum.specfit
        self.order = order
        fitp = np.zeros(self.order+1)
        self.spectofit = self.Spectrum.data+self.basespec
        self.OKmask = (self.spectofit==self.spectofit)
        if exclude == 'interactive' or interactive:
            self.excludemask[:] = True
            self.excludevelo = []
            self.excludepix  = []
            self.click = self.specplotter.axis.figure.canvas.mpl_connect('button_press_event',self.selectregion)
        else:
            if excludefit and specfit.modelpars is not None:
                #vlo = self.specplotter.specfit.modelpars[1] - 2*self.specplotter.specfit.modelpars[2]
                #vhi = self.specplotter.specfit.modelpars[1] + 2*self.specplotter.specfit.modelpars[2]
                #exclude = [np.argmin(abs(self.Spectrum.xarr-vlo)),argmin(abs(self.Spectrum.xarr-vhi))]
                specfit.fullsizemodel() # make sure the spectrum is the right size
                self.excludemask = abs(specfit.model) > exclusionlevel*abs(min(specfit.modelpars[0::3]))
            else:
                self.excludemask[:] = False
            self.dofit(exclude=exclude,annotate=annotate,**kwargs)
        if save: self.savefit()

    def dofit(self, exclude=None, excludeunits='velo', annotate=False,
            subtract=True, **kwargs):
        """
        Do the baseline fitting and save and plot the results.

        Can specify a region to exclude using velocity units or pixel units
        """
        if exclude is not None and excludeunits in ['velo','km/s','wavelength','frequency']:
            if len(exclude) % 2 == 0:
                self.excludevelo = exclude
                self.excludepix = []
                for vl,vu in zip(exclude[::2],exclude[1::2]):
                    xl = np.argmin(abs(self.Spectrum.xarr-vl))
                    xu = np.argmin(abs(self.Spectrum.xarr-vu))
                    if xl > xu: xl,xu=xu,xl
                    self.excludemask[xl:xu] = True
                    self.excludepix += [xl,xu]
        elif excludeunits in ['pix','pixel','chan','channel']:
            if len(exclude) % 2 == 0:
                self.excludepix = []
                for xl,xu in zip(exclude[::2],exclude[1::2]):
                    if xl > xu: xl,xu=xu,xl
                    self.excludemask[xl:xu] = True
                    self.excludepix += [xl,xu]
        self.basespec, self.baselinepars = self._baseline(
                self.spectofit,
                xarr=self.Spectrum.xarr,
                order=self.order, exclude=None, 
                mask=(True-self.OKmask)+self.excludemask,
                **kwargs)
        #self.basespec = np.poly1d(self.baselinepars)(self.Spectrum.xarr)
        if subtract:
            self.Spectrum.data -= self.basespec
            self.subtracted = True
        else:
            self.subtracted = False

        if self.specplotter.axis is not None:
            self.plot_baseline()

    def plot_baseline(self, annotate=True, plotcolor='orange'):
        if self.specplotter.axis is not None: 
            [self.specplotter.axis.lines.remove(p) for p in self.specplotter.axis.lines]
        if self.specplotter.errorplot is not None: 
            [self.specplotter.axis.collections.remove(p) for p in self.specplotter.errorplot if isinstance(p,matplotlib.collections.PolyCollection)]
            [self.specplotter.axis.lines.remove(p) for p in self.specplotter.errorplot if isinstance(p,matplotlib.lines.Line2D)]
        self.specplotter.plot(**self.specplotter.plotkwargs)
        if self.subtracted is False:
            self.specplotter.axis.plot(self.Spectrum.xarr,self.basespec,color=plotcolor)
        self.specplotter.axis.set_ylim(
                abs(self.Spectrum.data[self.OKmask].min())*1.1*np.sign(self.Spectrum.data[self.OKmask].min()),
                abs(self.Spectrum.data[self.OKmask].max())*1.1*np.sign(self.Spectrum.data[self.OKmask].max()))
        if annotate: self.annotate() # refreshes automatically
        elif self.specplotter.autorefresh: self.specplotter.refresh()

    def selectregion(self,event,annotate=False):
        """
        select regions for baseline fitting
        """
        if event.button == 1:
            if self.nclicks_b1 == 0:
                self.bx1 = np.argmin(abs(event.xdata-self.Spectrum.xarr))
                self.excludevelo += [self.Spectrum.xarr]
                self.excludepix  += [self.bx1]
                self.nclicks_b1 += 1
            elif self.nclicks_b1 == 1:
                self.bx2 = np.argmin(abs(event.xdata-self.Spectrum.xarr))
                self.nclicks_b1 -= 1
                if self.bx1 > self.bx2: self.bx1,self.bx2 = self.bx2,self.bx1
                self.fitregion += self.specplotter.axis.plot(
                        self.Spectrum.xarr[self.bx1:self.bx2],
                        self.Spectrum.data[self.bx1:self.bx2]+self.specplotter.offset,
                        drawstyle='steps-mid',
                        color='g',alpha=0.5)
                self.specplotter.refresh()
                self.excludemask[self.bx1:self.bx2] = False
                self.excludevelo += [self.Spectrum.xarr]
                self.excludepix  += [self.bx2]
        if event.button in [2,3]:
            disconnect(self.click)
            self.dofit(exclude=None,annotate=annotate)
            for p in self.fitregion:
                p.set_visible(False)
                self.specplotter.axis.lines.remove(p)
            self.fitregion=[] # I should be able to just remove from the list... but it breaks the loop...
            self.specplotter.refresh()

    def annotate(self,loc='upper left'):
        bltext = "bl: $y=$"+"".join(["$%+6.3gx^{%i}$" % (f,self.order-i)
            for i,f in enumerate(self.baselinepars)])
        #self.blleg = text(xloc,yloc     ,bltext,transform = self.specplotter.axis.transAxes)
        self.clearlegend()
        pl = matplotlib.collections.CircleCollection([0],edgecolors=['k'])
        self.blleg = self.specplotter.axis.legend(
                (pl,),
                (bltext,),loc=loc,markerscale=0.001,
                borderpad=0.1, handlelength=0.1, handletextpad=0.1
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

    def _baseline(self,spectrum,xarr=None,xmin='default',xmax='default',order=1,quiet=True,exclude=None,
            mask=None):
        """
        Subtract a baseline from a spectrum
        If xmin,xmax are not specified, defaults to ignoring first and last 10% of spectrum
        *unless* order > 1, in which case ignoring the ends tends to cause strange effects

        exclude is a set of start/end indices to ignore when baseline fitting
        (ignored by setting error to infinite in fitting procedure)
        """
        if xmin == 'default':
            if order <= 1: xmin = np.floor( spectrum.shape[-1]*0.1 )
            else:          xmin = 0
        elif xmin is None:
            xmin = 0
        if xmax == 'default':
            if order <= 1: xmax = np.ceil( spectrum.shape[-1]*0.9 )
            else:          xmax = spectrum.shape[-1]
        elif xmax is None:
            xmax = spectrum.shape[-1]
        
        pguess = [1]*(order+1)

        if xarr is None:
            xarr = np.indices(spectrum.shape).squeeze()

        subxarr = xarr[xmin:xmax]
        def mpfitfun(data,err):
            def f(p,fjac=None): return [0,np.ravel((np.poly1d(p)(subxarr)-data)/err)]
            return f

        err = np.ones(spectrum.shape)
        if exclude is not None:
            err[exclude[0]:exclude[1]] = 1e10
        if mask is not None:
            if mask.dtype.name != 'bool': mask = mask.astype('bool')
            err[mask] = 1e10
            spectrum[mask] = 0
        if (spectrum!=spectrum).sum() > 0:
            print "There is an error in baseline: some values are NaN"
            import pdb; pdb.set_trace()

        import mpfit
        mp = mpfit.mpfit(mpfitfun(spectrum[xmin:xmax],err[xmin:xmax]),xall=pguess,quiet=quiet)
        fitp = mp.params
        bestfit = np.poly1d(fitp)(xarr).squeeze()

        return bestfit,fitp
