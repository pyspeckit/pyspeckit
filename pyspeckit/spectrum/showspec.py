"""
showspec is my homegrown spectrum plotter, meant to somewhat follow STARLINK's
SPLAT and have functionality similar to GAIA, but with an emphasis on producing
publication-quality plots (which, while splat may do, it does unreproducibly)


TO DO:
    -add spectrum arithmetic tools
        (as is, you can use numpy.interp with sp.vind and sp.spectrum pretty
        easily)
    -implement other fitters
        -e.g., NH3 hyperfine, Voigt
    -move to object-oriented pylab/pyplot implementation (for bulk / speedup work)
    -allow for non-plotting fitting work (curious... I've never needed that yet)
    -Equivalent Width measurement without gaussian fitting
        -probably should be part of the baseline code
    -write documentation other people can read

"""



import math

import pylab
#from pylab import indices,figure,clf,savefig,plot,legend,text,axes,title,imshow,connect,get_current_fig_manager
from pylab import *
import matplotlib

from mpfit import mpfit

from collapse_gaussfit import *
from ratosexagesimal import *
import pyfits
import gaussfitter

from numpy import isnan
#from mad import MAD,nanmedian

def steppify(arr,isX=False,interval=0,sign=+1.0):
    """
    *support function*
    Converts an array to double-length for step plotting
    """
    if isX and interval==0:
        interval = abs(arr[1]-arr[0]) / 2.0
    newarr = array(zip(arr-sign*interval,arr+sign*interval)).ravel()
    return newarr

class SpecPlotter:
    """
    SpecPlotter class.  Takes in a spectrum or data cube, plotting properties,
    and a velocity axis determination function.  Look at splat_1d for a wrapper
    that might actually be useful.

    Whew, needs more documentation
    """

    def __init__(self,  cube,  axis=None,  xtol=None,  ytol=None, vconv=lambda
            x: x, xtora=lambda x: x, ytodec=lambda x: x, specname=None,
            dv=None, color='k', hdr=None, errspec=None,  maskspec=None,
            fig=None, fignum=1, clear=True, title=None, xunits='km/s',
            erralpha=0.2, ivconv=None, autorefresh=True, reffreq=None,
            gaiafignum=0, gaiafig=None, clickid=None, **kwargs ):
      self.vconv = vconv
      self.xtora = xtora
      self.ytodec = ytodec
      self.cube = cube # where(numpy.isnan(cube),0,cube)
      if len(self.cube.shape) > 1:
          self.spectrum = cube[:,0,0] # spectrum is what's plotted; cube is the "raw data"
      else:
          self.spectrum = cube # spectrum is what's plotted; cube is the "raw data"
      self.specname=specname
      self.dv=dv
      self.reffreq=reffreq
      self.scale=1.0
      self.units='K'
      self.xunits=xunits
      self.voff=0.0
      self.offset=0.0
      self.continuum=0.0
      self.errspec = errspec
      self.erralpha=erralpha
      self.plotcolor=color
      self.specfit = Specfit(self)
      self.fitspec = self.specfit
      self.baseline = Baseline(self)
      #self.fft = FFT(self)
      #self.psd = self.fft.psd
      self.vmin=None
      self.vmax=None
      self.title=title
      self.ivconv=ivconv
      self.autorefresh=autorefresh
      self.spectrumplot=None
      self.errorplot=None
      self.gaiafignum = gaiafignum
      self.gaiafig = gaiafig
      self.clickid = clickid
      self.plotkwargs = kwargs
      if maskspec is not None:
          self.maskspec = maskspec
      else:
          self.maskspec = zeros(self.cube.shape)
      self.linecollections =[]
      self.texts =[]
      if hdr: self.header = hdr
  
      # figure out where to put the plot
      if fig is None and axis is None:
          fig=figure(fignum)
          if clear: fig.clf()
          self.axis = fig.gca()
      elif fig is None and axis is None:
          self.axis = pylab.gca()
      elif fig is not None and axis is None:
          if clear: fig.clf()
          self.axis = fig.gca()
      elif fig is None and axis is not None:
          self.axis = axis
      else: # if figure and axis are both set, just use axis
          self.axis = axis
      if clear: self.axis.clear()
  
    def __call__(self, event):
      """
      Connects map cube to specplotter...
      """
      if event.inaxes:
        clickX = event.xdata
        clickY = event.ydata
        tb = get_current_fig_manager().toolbar
        #if ((self.axis is None) or (self.axis==event.inaxes)) and tb.mode=='':
        if event.button==1 and tb.mode=='':
            print "OverPlotting spectrum from point %i,%i" % (clickX,clickY)
            self.plotspec(clickY,clickX,button=event.button,cube=True)
        elif event.button==2:
            print "Plotting spectrum from point %i,%i" % (clickX,clickY)
            self.plotspec(clickY,clickX,button=event.button,cube=True,clear=True)
        elif event.button==3:
            print "Disconnecting GAIA-like tool"
            self.gaiafig.canvas.mpl_disconnect(self.clickid)
        else:
            print "Call failed for some reason: "
            print "event: ",event
  
    def plotspec(self, i=0, j=0, cube=False, title=None,
            button=1, clear=False,color=None, continuum=None,
            axis=None, offset=None, scale=None, voff=None, vmin=None,
            vmax=None, units=None, xunits=None, erralpha=None, plotpix=False,
            errstyle='fill', autorefresh=None, **kwargs):
      """
      Plot a spectrum
      Originally written to plot spectra from data cubes, hence the i,j parameter
      to specify the location in the cube
      """
  
      if kwargs.has_key('fignum'): kwargs.pop('fignum')  # HACK because I want __init__ to accept different kwargs
      if kwargs.has_key('fig'): kwargs.pop('fig')        # is there a better workaround?
      if scale  is not None: self.scale = scale
      if units  is not None: self.units = units
      if xunits is not None: self.xunits= xunits
      if voff   is not None: self.voff  = voff
      if offset is not None: self.offset= offset
      if continuum is not None: self.continuum= continuum
      if color  is not None: self.plotcolor=color
      if erralpha is not None: self.erralpha= erralpha
      if vmax  is not None: self.vmax = vmax
      if vmin  is not None: self.vmin = vmin
      if title is not None: self.title = title
      if autorefresh is not None: self.autorefresh = autorefresh
      if axis is None: axis=self.axis # allow spectrum to be plotted on other axis
      if clear: axis.clear()
  
      if plotpix:
          self.vind = arange(self.cube.shape[0])
      else:
          self.vind = self.vconv(arange(self.cube.shape[0])) + self.voff
  
      if kwargs.has_key('fignum'): kwargs.pop('fignum')
      if kwargs.has_key('linewidth'):
          linewidth = kwargs.pop('linewidth')
      else:
          linewidth=0.5
  
      if cube or len(self.cube.shape) == 3:
          self.spectrum = self.cube[:,i,j]*self.scale+self.continuum-self.baseline.basespec
          self.spectrumplot = axis.plot(self.vind,self.spectrum+self.offset,color=self.plotcolor,
                  linestyle='steps-mid',linewidth=linewidth,
                  **kwargs)
      else:
          if self.maskspec.sum() > 0:
              nanmask = where(self.maskspec,numpy.nan,1)
              self.spectrum = self.cube*self.scale*nanmask+self.continuum-self.baseline.basespec
              self.spectrumplot = axis.plot(self.vind,self.spectrum+self.offset,color=self.plotcolor,
                      linestyle='steps-mid',linewidth=linewidth,
                      **kwargs)
          else:
              self.spectrum = self.cube*self.scale+self.continuum-self.baseline.basespec
              self.spectrumplot = axis.plot(self.vind,self.spectrum+self.offset,color=self.plotcolor,
                      linestyle='steps-mid',linewidth=linewidth,
                      **kwargs)
          if self.errspec is not None:
              if errstyle == 'fill':
                  self.errorplot = [axis.fill_between(steppify(self.vind,isX=True,sign=sign(self.dv)),
                          steppify(self.spectrum+self.offset-self.errspec*self.scale),
                          steppify(self.spectrum+self.offset+self.errspec*self.scale),
                          facecolor=self.plotcolor, alpha=self.erralpha, **kwargs)]
              elif errstyle == 'bars':
                  self.errorplot = axis.errorbar(self.vind, self.spectrum+self.offset,
                          yerr=self.errspec*self.scale, ecolor=self.plotcolor, fmt=None,
                          **kwargs)
  
      if vmin is not None: xlo = self.vmin
      else: xlo=self.vind.min()
      if vmax is not None: xhi = self.vmax
      else: xhi=self.vind.max()
      axis.set_xlim(xlo,xhi)
  
      if self.title is not None:
          axis.set_title(self.title)
      elif self.xtora and self.ytodec:
          axis.set_title("Spectrum at %s %s" %
                  (ratos(self.xtora(i)),dectos(self.ytodec(j))))
      elif self.specname:
          axis.set_title("Spectrum of %s" % self.specname)
      if isinstance(self.xunits,str):
          axis.set_xlabel(self.xunits)
      else:
          axis.set_xlabel("V$_{LSR}$ (km s$^{-1}$)")
          self.xunits = 'km/s'
      if units in ['Ta*','Tastar','K']:
        axis.set_ylabel("$T_A^*$ (K)")
      elif units == 'mJy':
        axis.set_ylabel("$S_\\nu$ (mJy)")
      elif units == 'Jy':
        axis.set_ylabel("$S_\\nu$ (Jy)")
      else:
        axis.set_ylabel(self.units)
      if self.autorefresh: self.refresh()
  
    def save(self,fname,**kwargs):
      """
      Save the current spectrum (useful for saving baselined data)
      """
      newfile = pyfits.PrimaryHDU(data=self.cube,header=self.header)
      newfile.writeto(fname,**kwargs)

    def savefig(self,fname,bbox_inches='tight',**kwargs):
        """
        simple wrapper of maplotlib's savefig.  
        """
        self.axis.figure.savefig(fname,bbox_inches=bbox_inches,**kwargs)
  
    def showlines(self,linefreqs,linenames,ctype='freq',cunit='hz',yscale=0.8,vofflines=0.0,
            voffunit='km/s',**kwargs):
        """
        Overplot vertical lines and labels at the frequencies (or velocities) of each line
  
        yscale - fraction of maximum at which to label
        """

        self.clearlines() 

        if ctype != 'freq':
            print "Sorry, non-frequency units not implemented yet."
            return

        speedoflight=2.99792458e5
        if self.reffreq and self.xunits in ('km/s','m/s'):
            linefreqs = -(array(linefreqs)-self.reffreq)/self.reffreq * speedoflight
  
        if 'hz' in cunit or 'Hz' in cunit:
            linefreqs *= (1.0 + vofflines / speedoflight)
        else:
            linefreqs += vofflines
      
        ymax = (self.spectrum[self.spectrum==self.spectrum]).max() 
        for lf,ln in zip(linefreqs,linenames):
            if lf < self.vind.max() and lf > self.vind.min():
                self.linecollections.append(vlines(lf,0,ymax,**kwargs))
                self.texts.append(text(lf,ymax*yscale,ln,rotation='vertical',**kwargs))

        if self.autorefresh: self.refresh()
  
    def clearlines(self):
        if len(self.texts) > 0:
            for T in self.texts:
                if T in self.axis.texts:
                    self.axis.texts.remove(T)
        if len(self.linecollections) > 0:
            for LC in self.linecollections:
                if LC in self.axis.collections: 
                    self.axis.collections.remove(LC)
  
    def refresh(self):
        self.axis.figure.canvas.draw()

class FFT:
    def __init__(self,specplotter,fignum=3,axis=None, color='k'):
        self.specplotter=specplotter
        if axis is None:
            self.fignum=fignum
            self.figure=figure(self.fignum)
            self.axis=gca()
        else:
            self.axis=axis
            self.figure=self.axis.figure
            self.fignum=None
        #self.axis.clear()
        self.color=color
        self.fftplot=None
        self.setspec()
        self.setshift()
        self.clear()

    def __call__(self,psd=False,shift=True):
        self.setspec()
        if psd:
            self.psd(shift=shift)
        else:
            self.fft(shift=shift)
    
    def fft(self,shift=True,logplot=False,**kwargs):
        self.clear()
        self.setshift(shift)
        if logplot: self.axis.set_yscale('log')
        else: self.axis.set_yscale('linear')
        self.fftspec = fft(self.spectofft)
        self.realfft = self.fftspec.real
        self.imagfft = self.fftspec.imag
        self.fftplot = self.axis.plot(self.shiftfunc(self.realfft),
                drawstyle='steps-mid',color=self.color,**kwargs)
        self.refresh()

    def psd(self,logplot=True,shift=True,**kwargs):
        self.clear()
        if logplot: self.axis.set_yscale('log')
        else: self.axis.set_yscale('linear')
        self.setshift(shift)
        self.psdspec = fft(self.spectofft) * fft(self.spectofft[::-1])
        self.psdreal = abs(self.psdspec)
        self.fftplot = self.axis.plot(self.shiftfunc(self.psdreal),
                drawstyle='steps-mid',color=self.color,**kwargs)
        if logplot: self.axis.set_yscale('log')
        else: self.axis.set_yscale('linear')
        self.refresh()

    def setshift(self,shift=True):
        if shift: self.shiftfunc = fftshift
        else: self.shiftfunc = lambda x: x

    def setspec(self):
        self.spectofft = copy(self.specplotter.spectrum)
        OKmask = (self.spectofft==self.spectofft)
        self.spectofft[(True-OKmask)] = 0

    def clear(self):
        if self.fftplot is not None:
            for p in self.fftplot:
                p.set_visible(False)
                if p in self.axis.lines: self.axis.lines.remove(p)
        self.axis.clear()
        self.refresh()

    def refresh(self):
        self.axis.figure.canvas.draw()

class PSD(FFT):
    def __call__(self,shift=True):
        self.setspec()
        self.setshift(shift)
        self.clear()
        self.psd()
        self.refresh()

class Baseline:
    def __init__(self,specplotter):
        self.baselinepars  = None
        self.order = None
        self.basespec = zeros(specplotter.spectrum.shape[0])
        self.excludemask = zeros(specplotter.spectrum.shape[0],dtype='bool')
        self.OKmask = ones(specplotter.spectrum.shape[0],dtype='bool')
        self.specplotter = specplotter
        self.blleg = None
        self.click = 0
        self.nclicks_b1 = 0
        self.nclicks_b2 = 0
        self.fitregion=[]
        self.excludevelo = []
        self.excludepix  = []

    def __call__(self, order=1, annotate=False, excludefit=False, save=True,
            exclude=None, exclusionlevel=0.01,
            interactive=False, **kwargs):
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
        specfit = self.specplotter.specfit
        self.order = order
        fitp = zeros(self.order+1)
        self.spectofit = self.specplotter.spectrum+self.basespec
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
                #exclude = [argmin(abs(self.specplotter.vind-vlo)),argmin(abs(self.specplotter.vind-vhi))]
                specfit.fullsizemodel() # make sure the spectrum is the right size
                self.excludemask = abs(specfit.model) > exclusionlevel*abs(min(specfit.modelpars[0::3]))
            else:
                self.excludemask[:] = False
            self.dofit(exclude=exclude,annotate=annotate,**kwargs)
        if save: self.savefit()

    def dofit(self, exclude=None, excludeunits='velo', annotate=False,
            **kwargs):
        """
        Do the baseline fitting and save and plot the results.

        Can specify a region to exclude using velocity units or pixel units
        """
        if exclude is not None and excludeunits in ['velo','km/s']:
            if len(exclude) % 2 == 0:
                self.excludevelo = exclude
                self.excludepix = []
                for vl,vu in zip(exclude[::2],exclude[1::2]):
                    xl = argmin(abs(self.specplotter.vind-vl))
                    xu = argmin(abs(self.specplotter.vind-vu))
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
        self.specplotter.spectrum, self.baselinepars = baseline(
                self.spectofit,
                xarr=self.specplotter.vind,
                order=self.order, exclude=None, 
                mask=(True-self.OKmask)+self.excludemask,
                **kwargs)
        self.basespec = poly1d(self.baselinepars)(self.specplotter.vind)
        if self.specplotter.spectrumplot is not None: 
            [self.specplotter.axis.lines.remove(p) for p in self.specplotter.spectrumplot]
        if self.specplotter.errorplot is not None: 
            [self.specplotter.axis.collections.remove(p) for p in self.specplotter.errorplot if isinstance(p,matplotlib.collections.PolyCollection)]
            [self.specplotter.axis.lines.remove(p) for p in self.specplotter.errorplot if isinstance(p,matplotlib.lines.Line2D)]
        self.specplotter.plotspec(**self.specplotter.plotkwargs)
        self.specplotter.axis.set_ylim(
                abs(self.specplotter.spectrum[self.OKmask].min())*1.1*sign(self.specplotter.spectrum[self.OKmask].min()),
                abs(self.specplotter.spectrum[self.OKmask].max())*1.1*sign(self.specplotter.spectrum[self.OKmask].max()))
        if annotate: self.annotate() # refreshes automatically
        elif self.specplotter.autorefresh: self.specplotter.refresh()

    def selectregion(self,event,annotate=False):
        """
        select regions for baseline fitting
        """
        if event.button == 1:
            if self.nclicks_b1 == 0:
                self.bx1 = argmin(abs(event.xdata-self.specplotter.vind))
                self.excludevelo += [self.specplotter.vind]
                self.excludepix  += [self.bx1]
                self.nclicks_b1 += 1
            elif self.nclicks_b1 == 1:
                self.bx2 = argmin(abs(event.xdata-self.specplotter.vind))
                self.nclicks_b1 -= 1
                if self.bx1 > self.bx2: self.bx1,self.bx2 = self.bx2,self.bx1
                self.fitregion += self.specplotter.axis.plot(
                        self.specplotter.vind[self.bx1:self.bx2],
                        self.specplotter.spectrum[self.bx1:self.bx2]+self.specplotter.offset,
                        drawstyle='steps-mid',
                        color='g',alpha=0.5)
                self.specplotter.refresh()
                self.excludemask[self.bx1:self.bx2] = False
                self.excludevelo += [self.specplotter.vind]
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
        if self.baselinepars is not None:
            for ii,p in enumerate(self.baselinepars):
                self.specplotter.header.update('BLCOEF%0.2i' % (ii),p,comment="Baseline power-law best-fit coefficient x^%i" % (self.order-ii-1))


class Specfit:

    def __init__(self,specplotter):
        self.model = None
        self.modelpars = None
        self.modelerrs = None
        self.modelplot = None
        self.guessplot = []
        self.fitregion = []
        self.ngauss = 0
        self.nclicks_b1 = 0
        self.nclicks_b2 = 0
        self.gx1 = 0
        self.gx2 = specplotter.spectrum.shape[0]
        self.guesses = []
        self.click = 0
        self.fitkwargs = {}
        self.auto = False
        self.autoannotate = True
        self.specplotter = specplotter
        self.gaussleg=None
        self.residuals=None
        self.setfitspec()
        #self.seterrspec()

    def __call__(self, interactive=False, usemoments=True, fitcolor='r',
            multifit=False, guesses=None, annotate=True, save=True,
            **kwargs):
        """
        Fit gaussians to a spectrum

        guesses = [height,amplitude,center,width]
        """
  
        self.fitcolor = fitcolor
        self.clear()

        self.ngauss = 0
        self.fitkwargs = kwargs
        if interactive:
            print "Left-click twice to select a fitting range, then middle-click twice to select a peak and width"
            self.nclicks_b1 = 0
            self.nclicks_b2 = 0
            self.guesses = []
            self.click = self.specplotter.axis.figure.canvas.mpl_connect('button_press_event',self.makeguess)
            self.autoannotate = annotate
        elif multifit:
            if guesses is None:
                print "You must input guesses when using multifit.  Also, baseline first!"
            else:
                self.guesses = guesses
                self.multifit()
                self.autoannotate = annotate
        else:
            #print "Non-interactive, 1D fit with automatic guessing"
            if self.specplotter.baseline.order is None:
                self.specplotter.baseline.order=0
                self.onedfit(usemoments=usemoments,annotate=annotate,**kwargs)
            else:
                self.onedfit(usemoments=usemoments,annotate=annotate,
                        vheight=False,height=0.0,**kwargs)
            if self.specplotter.autorefresh: self.specplotter.refresh()
        if save: self.savefit()

    def seterrspec(self,usestd=None,useresiduals=True):
        if self.specplotter.errspec is not None and not usestd:
            self.errspec = self.specplotter.errspec
        elif self.residuals is not None and useresiduals: 
            self.errspec = ones(self.spectofit.shape[0]) * self.residuals.std()
        else: self.errspec = ones(self.spectofit.shape[0]) * self.spectofit.std()

    def setfitspec(self):
        self.spectofit = copy(self.specplotter.spectrum)
        OKmask = (self.spectofit==self.spectofit)
        self.spectofit[(True-OKmask)] = 0
        self.seterrspec()
        self.errspec[(True-OKmask)] = 1e10

    def multifit(self):
        self.ngauss = len(self.guesses)/3
        self.setfitspec()
        if self.fitkwargs.has_key('negamp'): self.fitkwargs.pop('negamp')
        mpp,model,mpperr,chi2 = gaussfitter.multigaussfit(
                self.specplotter.vind[self.gx1:self.gx2], 
                self.spectofit[self.gx1:self.gx2], 
                err=self.errspec[self.gx1:self.gx2],
                ngauss=self.ngauss,
                params=self.guesses,
                **self.fitkwargs)
        self.chi2 = chi2
        self.dof  = self.gx2-self.gx1-self.ngauss*3
        self.model = model
        self.modelpars = mpp.tolist()
        self.modelerrs = mpperr.tolist()
        self.modelplot = self.specplotter.axis.plot(
                self.specplotter.vind[self.gx1:self.gx2],
                self.model+self.specplotter.offset, color=self.fitcolor, linewidth=0.5)
        self.residuals = self.spectofit[self.gx1:self.gx2] - self.model
        if self.autoannotate:
            self.annotate()
    
    def onedfit(self, usemoments=True, annotate=True, vheight=True, height=0, negamp=None,**kwargs):
        self.ngauss = 1
        self.auto = True
        self.setfitspec()
        if usemoments: # this can be done within gaussfit but I want to save them
            self.guesses = gaussfitter.onedmoments(
                    self.specplotter.vind[self.gx1:self.gx2],
                    self.spectofit[self.gx1:self.gx2],
                    vheight=vheight,negamp=negamp,**kwargs)
            if vheight is False: self.guesses = [height]+self.guesses
        else:
            if negamp: self.guesses = [height,-1,0,1]
            else:  self.guesses = [height,1,0,1]
        mpp,model,mpperr,chi2 = gaussfitter.onedgaussfit(
                self.specplotter.vind[self.gx1:self.gx2],
                self.spectofit[self.gx1:self.gx2],
                err=self.errspec[self.gx1:self.gx2],
                vheight=vheight,
                params=self.guesses,
                **self.fitkwargs)
        self.chi2 = chi2
        self.dof  = self.gx2-self.gx1-self.ngauss*3-vheight
        if vheight: 
            self.specplotter.baseline.baselinepars = mpp[:1] # first item in list form
            self.model = model - mpp[0]
        else: self.model = model
        self.residuals = self.spectofit[self.gx1:self.gx2] - self.model
        self.modelpars = mpp[1:].tolist()
        self.modelerrs = mpperr[1:].tolist()
        self.modelplot = self.specplotter.axis.plot(
                self.specplotter.vind[self.gx1:self.gx2],
                self.model+self.specplotter.offset, color=self.fitcolor, linewidth=0.5)
        if annotate:
            self.annotate()
            if vheight: self.specplotter.baseline.annotate()

    def fullsizemodel(self):
        """
        If the gaussian was fit to a sub-region of the spectrum,
        expand it (with zeros) to fill the spectrum.  You can 
        always recover the original by:
        origmodel = model[gx1:gx2]
        """

        if self.model.shape != self.specplotter.spectrum.shape:
            temp = zeros(self.specplotter.spectrum.shape)
            temp[self.gx1:self.gx2] = self.model
            self.model = temp
            self.residuals = self.spectofit - self.model

    def plotresiduals(self,fig=None,axis=None,clear=True,**kwargs):
        """
        Plot residuals of the fit.  Specify a figure or
        axis; defaults to figure(2).

        kwargs are passed to matplotlib plot
        """
        if axis is None:
            fig=figure(2)
            self.residualaxis = gca()
            if clear: self.residualaxis.clear()
        else:
            self.residualaxis = axis
            if clear: self.residualaxis.clear()
        self.residualplot = self.residualaxis.plot(self.specplotter.vind[self.gx1:self.gx2],
                self.residuals,drawstyle='steps-mid',
                linewidth=0.5, color='k', **kwargs)
        if self.specplotter.vmin is not None and self.specplotter.vmax is not None:
            self.residualaxis.set_xlim(self.specplotter.vmin,self.specplotter.vmax)
        self.residualaxis.figure.canvas.draw()

    def annotate(self,loc='upper right'):
        #text(xloc,yloc     ,"c=%g" % self.modelpars[1],transform = self.specplotter.axis.transAxes)
        #text(xloc,yloc-0.05,"w=%g" % self.modelpars[2],transform = self.specplotter.axis.transAxes)
        #text(xloc,yloc-0.10,"a=%g" % self.modelpars[0],transform = self.specplotter.axis.transAxes)
        self.clearlegend()
        pl = matplotlib.collections.CircleCollection([0],edgecolors=['k'])
        self.gaussleg = self.specplotter.axis.legend(
                tuple([pl]*3*self.ngauss),
                tuple(flatten(
                    [("c%i=%6.4g $\\pm$ %6.4g" % (jj,self.modelpars[1+jj*3],self.modelerrs[1+jj*3]),
                      "w%i=%6.4g $\\pm$ %6.4g" % (jj,self.modelpars[2+jj*3],self.modelerrs[2+jj*3]),
                      "a%i=%6.4g $\\pm$ %6.4g" % (jj,self.modelpars[0+jj*3],self.modelerrs[0+jj*3]))
                      for jj in range(self.ngauss)])),
                loc=loc,markerscale=0.01,
                borderpad=0.1, handlelength=0.1, handletextpad=0.1
                )
        self.gaussleg.draggable(True)
        self.specplotter.axis.add_artist(self.gaussleg)
        if self.specplotter.autorefresh: self.specplotter.refresh()

    def selectregion(self,event):
        if self.nclicks_b1 == 0:
            self.gx1 = argmin(abs(event.xdata-self.specplotter.vind))
            self.nclicks_b1 += 1
        elif self.nclicks_b1 == 1:
            self.gx2 = argmin(abs(event.xdata-self.specplotter.vind))
            self.nclicks_b1 -= 1
            if self.gx1 > self.gx2: self.gx1,self.gx2 = self.gx2,self.gx1
            if abs(self.gx1-self.gx2) > 3: # can't fit w/ fewer data than pars
                self.fitregion = self.specplotter.axis.plot(
                        self.specplotter.vind[self.gx1:self.gx2],
                        self.specplotter.spectrum[self.gx1:self.gx2]+self.specplotter.offset,
                        drawstyle='steps-mid',
                        color='c')
                if self.guesses == []:
                    self.guesses = gaussfitter.onedmoments(
                            self.specplotter.vind[self.gx1:self.gx2],
                            self.spectofit[self.gx1:self.gx2],
                            vheight=0)
                    self.ngauss = 1
                    self.auto = True
            else:
                print "Fitting region is too small (channels %i:%i).  Try again." % (self.gx1,self.gx2)

    def guesspeakwidth(self,event):
        if self.nclicks_b2 % 2 == 0:
            if self.auto:
                self.guesses[:2] = [event.ydata,event.xdata]
            else:
                self.guesses += [event.ydata,event.xdata,1]
                self.ngauss += 1
            self.nclicks_b2 += 1
            self.guessplot += [self.specplotter.axis.scatter(event.xdata,event.ydata,marker='x',c='r')]
        elif self.nclicks_b2 % 2 == 1:
            self.guesses[-1] = abs(event.xdata-self.guesses[-2])
            self.nclicks_b2 += 1
            self.guessplot += self.specplotter.axis.plot([event.xdata,
                2*self.guesses[-2]-event.xdata],[event.ydata]*2,
                color='r')
            if self.auto:
                self.auto = False
            if self.nclicks_b2 / 2 > self.ngauss:
                print "There have been %i middle-clicks but there are only %i gaussians" % (self.nclicks_b2,self.ngauss)
                self.ngauss += 1

    def clear(self,legend=True):
        if self.modelplot is not None:
            for p in self.modelplot:
                p.set_visible(False)
        if legend: self.clearlegend()

    def makeguess(self,event):
        if event.button == 1:
            self.selectregion(event)
        elif event.button == 2:
            self.guesspeakwidth(event)
        elif event.button == 3:
            disconnect(self.click)
            if self.ngauss > 0:
                print len(self.guesses)/3," Guesses: ",self.guesses," X channel range: ",self.gx1,self.gx2
                if len(self.guesses) % 3 == 0:
                    self.multifit()
                    for p in self.guessplot + self.fitregion:
                        p.set_visible(False)
                else: 
                    print "error, wrong # of pars"
        if self.specplotter.autorefresh: self.specplotter.refresh()

    def clearlegend(self):
        if self.gaussleg is not None: 
            self.gaussleg.set_visible(False)
            if self.gaussleg in self.specplotter.axis.artists:
                self.specplotter.axis.artists.remove(self.gaussleg)
        if self.specplotter.autorefresh: self.specplotter.refresh()
    
    def savefit(self):
        if self.modelpars is not None:
            for ii,p in enumerate(self.modelpars):
                if ii % 3 == 0: self.specplotter.header.update('AMP%1i' % (ii/3),p,comment="Gaussian best fit amplitude #%i" % (ii/3))
                if ii % 3 == 1: self.specplotter.header.update('CEN%1i' % (ii/3),p,comment="Gaussian best fit center #%i" % (ii/3))
                if ii % 3 == 2: self.specplotter.header.update('WID%1i' % (ii/3),p,comment="Gaussian best fit width #%i" % (ii/3))


def mapplot(plane,cube,vconv=lambda x: x,xtora=lambda x: x,ytodec=lambda x: x, gaiafignum=0, specfignum=1):
    
    gaiafig = figure(gaiafignum)
    gaiafig.clf()
    gaiaax = gaiafig.add_subplot(111)
    gaiaax.imshow(plane)

    sp = SpecPlotter(cube, vconv=vconv, xtora=xtora, ytodec=ytodec,
            gaiafignum=gaiafignum, fignum=specfignum, gaiafig=gaiafig)
    sp.clickid = gaiafig.canvas.mpl_connect('button_press_event',sp)
    #connect('button_press_event',sp)

 
def splat_3d(filename,xi=0,yi=0,vmin=None,vmax=None,button=1,dobaseline=False,exclude=None,
        smooth=None,smoothto=None,smoothtype='gaussian',order=1,savepre=None,**kwargs):
    """
    Inputs:
        vmin,vmax - range over which to baseline and plottransform = ax.transAxes
        exclude - (internal) range to exclude from baseline fit
    """
    dv,v0,p3,hdr,cube,xtora,ytodec,vconv,xunits,conversion_factor,units = open_3d(filename)

    if units is None: units="UNITS"
    if xunits is None: xunits="km/s"
    if conversion_factor == 0 or conversion_factor is None: conversion_factor=1.0
    sp = splat_1d(vpars=[dv, v0, p3], hdr=hdr, spec=cube[:, yi, xi],
            xtora=xtora, ytodec=ytodec, vconv=vconv, units=units,
            conversion_factor=conversion_factor, xunits=xunits, **kwargs)

    sp.cube = cube

    return sp

def gaia(filename,estimator='max',axis=0):
    f = pyfits.open(filename)
    hdr = f[0].header
    cube = f[0].data
    dv,v0,p3 = hdr.get('CD3_3'),hdr.get('CRVAL3'),hdr.get('CRPIX3')
    dr,r0,p1 = hdr.get('CD1_1'),hdr.get('CRVAL1'),hdr.get('CRPIX1')
    dd,d0,p2 = hdr.get('CD2_2'),hdr.get('CRVAL2'),hdr.get('CRPIX2')
    if dv is None: dv = hdr.get('CDELT3')
    if dr is None: dr = hdr.get('CDELT1')
    if dd is None: dd = hdr.get('CDELT2')
    xtora = lambda x: (x-p1+1)*dr+r0    # convert pixel coordinates to RA/Dec/Velocity
    ytodec = lambda y: (y-p2+1)*dd+d0
    vconv = lambda v: (v-p3+1)*dv+v0

    if axis > 0:
        cube = cube.swapaxes(0,axis)

    if estimator == 'max':
        p = where(isnan(cube),0,cube).max(axis=0)
    elif estimator == 'int':
        p = where(isnan(cube),0,cube).sum(axis=0) * dv
    elif estimator == 'intdivmax':
        cut = MAD(cube.ravel()) + nanmedian(cube.ravel())
        if cut < 0:
            cut = 0
        m = where(isnan(cube),0,cube).max(axis=0)
        i = where(isnan(cube),0,cube).sum(axis=0) * dv
        p = where(i<0,0,i)/where(m<=cut,numpy.inf,m)
    elif estimator[-5:] == ".fits":
        p = pyfits.open(estimator)[0].data

    mapplot(p,cube,vconv,xtora,ytodec)

def baseline_file(filename,outfilename,vmin=None,vmax=None,order=1,crop=False):
    f = pyfits.open(filename)
    hdr = f[0].header
    cube = f[0].data.squeeze()
    dv,v0,p3 = hdr.get('CD3_3'),hdr.get('CRVAL3'),hdr.get('CRPIX3')
    dr,r0,p1 = hdr.get('CD1_1'),hdr.get('CRVAL1'),hdr.get('CRPIX1')
    dd,d0,p2 = hdr.get('CD2_2'),hdr.get('CRVAL2'),hdr.get('CRPIX2')
    if dv is None: dv = hdr.get('CDELT3')
    if dr is None: dr = hdr.get('CDELT1')
    if dd is None: dd = hdr.get('CDELT2')
    vconv = lambda v: (v-p3+1)*dv+v0
    varr = vconv(arange(cube.shape[-1]))
    if vmin is None: argvmin = None
    else: argvmin = argmin(abs(varr-vmin))
    if vmax is None: argvmax = None
    else: argvmax = argmin(abs(varr-vmax))

    bspec,bfit = baseline(cube,vmin=argvmin,vmax=argvmax,order=order)


def baseline(spectrum,xarr=None,xmin='default',xmax='default',order=1,quiet=True,exclude=None,
        mask=None):
    """
    Subtract a baseline from a spectrum
    If xmin,xmax are not specified, defaults to ignoring first and last 10% of spectrum
    *unless* order > 1, in which case ignoring the ends tends to cause strange effects

    exclude is a set of start/end indices to ignore when baseline fitting
    (ignored by setting error to infinite in fitting procedure)
    """
    if xmin == 'default':
        if order <= 1: xmin = floor( spectrum.shape[-1]*0.1 )
        else:          xmin = 0
    elif xmin is None:
        xmin = 0
    if xmax == 'default':
        if order <= 1: xmax = ceil( spectrum.shape[-1]*0.9 )
        else:          xmax = spectrum.shape[-1]
    elif xmax is None:
        xmax = spectrum.shape[-1]
    
    pguess = [1]*(order+1)

    if xarr is None:
        xarr = indices(spectrum.shape).squeeze()

    subxarr = xarr[xmin:xmax]
    def mpfitfun(data,err):
        def f(p,fjac=None): return [0,numpy.ravel((poly1d(p)(subxarr)-data)/err)]
        return f

    err = ones(spectrum.shape)
    if exclude is not None:
        err[exclude[0]:exclude[1]] = 1e10
    if mask is not None:
        if mask.dtype.name != 'bool': mask = mask.astype('bool')
        err[mask] = 1e10
        spectrum[mask] = 0
    if (spectrum!=spectrum).sum() > 0:
        print "There is an error in baseline: some values are NaN"
        import pdb; pdb.set_trace()

    mp = mpfit(mpfitfun(spectrum[xmin:xmax],err[xmin:xmax]),xall=pguess,quiet=quiet)
    fitp = mp.params
    bestfit = poly1d(fitp)(xarr).squeeze()

    return (spectrum-bestfit),fitp

def open_3d(filename):
    f = pyfits.open(filename)
    hdr = f[0].header
    cube = f[0].data
    if len(cube.shape) == 4: cube=cube[0,:,:,:]
    #cube = reshape(cube.mean(axis=2).mean(axis=1),[cube.shape[0],1,1])
    dv,v0,p3 = hdr.get('CD3_3'),hdr.get('CRVAL3'),hdr.get('CRPIX3')
    dr,r0,p1 = hdr.get('CD1_1'),hdr.get('CRVAL1'),hdr.get('CRPIX1')
    dd,d0,p2 = hdr.get('CD2_2'),hdr.get('CRVAL2'),hdr.get('CRPIX2')
    if dv is None: dv = hdr.get('CDELT3')
    if dr is None: dr = hdr.get('CDELT1')
    if dd is None: dd = hdr.get('CDELT2')
    xtora = lambda x: (x-p1+1)*dr+r0    # convert pixel coordinates to RA/Dec/Velocity
    ytodec = lambda y: (y-p2+1)*dd+d0
    vconv = lambda v: (v-p3+1)*dv+v0
    if hdr.get('CUNIT3') in ['m/s','M/S']:
        conversion_factor = 1000.0
        xunits = 'km/s' # change to km/s because you're converting units
    else:
        xunits = hdr.get('CUNIT3')
        if xunits in ("hz","Hz"):
            print "Converting from Hz to GHz"
            xunits = "GHz"
            conversion_factor = 1.0e9
        else:
            conversion_factor = 1.0
    units = hdr.get('BUNIT')

    return dv,v0,p3,hdr,cube,xtora,ytodec,vconv,xunits,conversion_factor,units

def open_1d(filename,specnum=0,wcstype='',errspecnum=None,maskspecnum=None):
    """
    Grabs all the relevant pieces of a 1d spectrum for plotting
    wcstype is the suffix on the WCS type to get to velocity/frequency/whatever
    """
    f = pyfits.open(filename)
    hdr = f[0].header
    spec = f[0].data
    errspec  = None
    maskspec = None
    if hdr.get('NAXIS') == 2:
        if errspecnum is not None:
            errspec = spec[errspecnum,:]
        if maskspecnum is not None:
            maskspec = spec[maskspecnum,:]
        if isinstance(specnum,list):
            spec = spec[specnum,:].mean(axis=0)
        elif isinstance(specnum,int):
            spec = spec[specnum,:]
        else:
            raise TypeError("Specnum is of wrong type (not a list of integers or an integer).  Type: %s" %
                    str(type(specnum)))
    elif hdr.get('NAXIS') > 2:
        raise ValueError("Too many axes for open_1d (splat_1d) - use cube instead")
    if hdr.get('CD1_1'+wcstype):
        dv,v0,p3 = hdr['CD1_1'+wcstype],hdr['CRVAL1'+wcstype],hdr['CRPIX1'+wcstype]
    else:
        dv,v0,p3 = hdr['CDELT1'+wcstype],hdr['CRVAL1'+wcstype],hdr['CRPIX1'+wcstype]
    if hdr.get('OBJECT'):
        specname = hdr.get('OBJECT')
    elif hdr.get('GLON') and hdr.get('GLAT'):
        specname = "%s %s" % (hdr.get('GLON'),hdr.get('GLAT'))
    else:
        specname = filename.rstrip(".fits")
    if hdr.get('CUNIT1'+wcstype) in ['m/s','M/S']:
        conversion_factor = 1000.0
        xunits = 'km/s' # change to km/s because you're converting units
    else:
        xunits = hdr.get('CUNIT1'+wcstype)
        if xunits in ("hz","Hz"):
            print "Converting from Hz to GHz"
            xunits = "GHz"
            conversion_factor = 1.0e9
        else:
            conversion_factor = 1.0
    vconv = lambda v: ((v-p3+1)*dv+v0)/conversion_factor
    xtora=None
    ytodec=None
    units = hdr.get('BUNIT').strip()
    if hdr.get('CTYPE1'+wcstype):
        xtype = hdr.get('CTYPE1'+wcstype)
    else:
        xtype = 'VLSR'
    if hdr.get('REFFREQ'+wcstype):
        reffreq = hdr.get('REFFREQ'+wcstype)
    else:
        reffreq = None

    return dv,v0,p3,conversion_factor,hdr,spec,vconv,xtora,ytodec,specname,units,xunits,errspec,maskspec,reffreq

def splat_1d(filename=None,vmin=None,vmax=None,button=1,dobaseline=False,
        exclude=None,smooth=None,order=1,savepre=None,vcrop=True,
        vconv=None,vpars=None,hdr=None,spec=None,xtora=None,ytodec=None,
        specname=None,quiet=True,specnum=0,errspecnum=None,wcstype='',
        offset=0.0, continuum=0.0, annotatebaseline=False, plotspectrum=True,
        smoothto=None, xunits=None, units=None, conversion_factor=None,
        smoothtype='gaussian',convmode='valid',maskspecnum=None,**kwargs):
    """
    Wrapper for specplotter creation.  Works nicely with 1D spectra with well-defined
    FITS headers (i.e., CRVAL1, CRPIX1, CDELT1, and optionally CUNIT1 and CTYPE1)

    This documentation needs to be updated a lot... I implemented a lot of features
    without documenting them, which was a mistake

    Inputs:
        vmin,vmax - range over which to baseline and plot
        exclude - (internal) range to exclude from baseline fit
        vcrop - will vmin/vmax crop out data, or just set the plot limits?
    """
    if (vpars and vconv and hdr and spec is not None and xtora and ytodec 
            and units and xunits and conversion_factor):
        dv,v0,p3 = vpars
        errspec = None
        maskspec = None
        reffreq = None
        if units is None and kwargs.has_key('units'): units = kwargs.pop('units')
    else:
        dv,v0,p3,conversion_factor,hdr,spec,vconv,xtora,ytodec,specname_file,units,xunits,errspec,maskspec,reffreq = \
                open_1d(filename,specnum=specnum,wcstype=wcstype,errspecnum=errspecnum,maskspecnum=maskspecnum)
        if specname is None: specname=specname_file
        if units is None and kwargs.has_key('units'): units = kwargs.pop('units')

    if type(continuum)==type('str'):
        if hdr.get(continuum) is not None:
            continuum = hdr.get(continuum)
        else:
            raise ValueError("Continuum specified but none present.")

    varr = vconv(arange(spec.shape[0]))
    if vmin is None or vcrop==False: argvmin = 0
    else: 
      argvmin = argmin(abs(varr-vmin))
      if dv > 0:
        hdr.update('CRPIX1'+wcstype,p3-argvmin)
    if vmax is None or vcrop==False: argvmax = spec.shape[0]
    else: 
      argvmax = argmin(abs(varr-vmax))
      if dv < 0:
        hdr.update('CRPIX1'+wcstype,p3-argvmax)

    if argvmin > argvmax:
        argvmin,argvmax = argvmax,argvmin
        #if exclude is not None: exclude = exclude[::-1]
    elif argvmin == argvmax:
        raise Exception("Error: no data in velocity range %g:%g for source %s."
                % (vmin,vmax,filename))

    # these lines were meant to automatically put "exclude" into velocity
    # units; this is now done in the baseline code
    #if exclude is not None:
    #    exclude[0] = argmin(abs(varr-exclude[0]))
    #    exclude[1] = argmin(abs(varr-exclude[1]))
    #    exclude = array(exclude) - argvmin

    vconv = lambda v: ((v-p3+argvmin+1)*dv+v0) / conversion_factor
    ivconv = lambda V: p3-1-argvmin+(V*conversion_factor-v0)/dv

    specplot = spec[argvmin:argvmax]
    if errspec is not None: errspec=errspec[argvmin:argvmax]
    if maskspec is not None: maskspec=maskspec[argvmin:argvmax]

    if smoothto:
        smooth = abs(smoothto/dv)

    if smooth:
        roundsmooth = round(smooth) # can only downsample by integers
        # change fitter first
        if smoothtype == 'hanning': 
            specplot = convolve(specplot,hanning(2+roundsmooth)/hanning(2+roundsmooth).sum(),convmode)[::roundsmooth]
            kernsize = smooth
            ones_sameshape = zeros(smooth+2)
            ones_sameshape[1:-1] = 1
        elif smoothtype == 'boxcar':
            specplot = convolve(specplot,ones(roundsmooth)/float(roundsmooth),convmode)[::roundsmooth]
            kernsize = roundsmooth
            ones_sameshape = ones(roundsmooth)
        elif smoothtype == 'gaussian':
            speclen = specplot.shape[0]
            xkern  = linspace(-1*smooth,smooth,smooth*3)
            kernel = exp(-xkern**2/(2*(smooth/sqrt(8*log(2)))**2))
            kernel /= kernel.sum()
            kernsize = len(kernel)
            specplot = convolve(specplot,kernel,convmode)[::roundsmooth] 
            ones_sameshape = zeros(roundsmooth*3)
            ones_sameshape[roundsmooth:-roundsmooth] = 1
        if errspec is not None: 
            errspec = sqrt(convolve(errspec**2,ones_sameshape,convmode)[::roundsmooth]) / float(roundsmooth)
        if maskspec is not None: 
            maskspec = array(convolve(maskspec,ones_sameshape,convmode)[::roundsmooth],dtype='bool')
            if maskspec.shape != specplot.shape: import pdb; pdb.set_trace()
        # this bit of code may also make sense, but I'm shifting the center pixel instead
        # b/c it's easier (?) to deal with velocity range
        #v0 += (abs(dv)*smooth - abs(dv))/2.0 # pixel center moves by half the original pixel size
        dv *= roundsmooth
        if convmode == 'same':
          newrefpix = (p3-argvmin)/roundsmooth  
        elif convmode == 'full':
          newrefpix = (p3-0.5-argvmin+kernsize/2.0)/roundsmooth  
        elif convmode == 'valid':
          newrefpix = (p3-0.5-argvmin-kernsize/2.0)/roundsmooth  
        # this was resolved by advanced guess-and check
        # but also, sort of makes sense: FITS refers to the *center* of a pixel.  You want to 
        # shift 1/2 pixel to the right so that the first pixel goes from 0 to 1
        vconv = lambda v: ((v-newrefpix)*dv+v0)/conversion_factor
        ivconv = lambda V: newrefpix+(V*conversion_factor-v0)/dv
        hdr.update('CRPIX1'+wcstype,newrefpix+1)
        hdr.update('CDELT1'+wcstype,dv)

    sp = SpecPlotter(specplot, vconv=vconv, xtora=xtora, ytodec=ytodec,
            specname=specname, dv=dv/conversion_factor, hdr=hdr, reffreq=reffreq,
            errspec=errspec, maskspec=maskspec, xunits=xunits, **kwargs)

    if plotspectrum:
        sp.plotspec(button=button, cube=False, vmin=vmin, vmax=vmax,
                units=units, offset=offset, continuum=continuum,
                **kwargs)

    if dobaseline: 
        sp.baseline(exclude=exclude,order=order,quiet=quiet,annotate=annotatebaseline)
    if plotspectrum: sp.refresh()
    
    if hdr.get('GLON') and hdr.get('GLAT'):
        sp.glon = hdr.get('GLON')
        sp.glat = hdr.get('GLAT')

    if savepre is not None:
        glon,glat = sp.glon,sp.glat
        if glat < 0: pm="" 
        else: pm = "+"
        savename = savepre + "G%07.3f%0s%07.3f_" % (glon,pm,glat) + hdr.get('MOLECULE').replace(' ','') + hdr.get('TRANSITI').replace(' ','')
        savefig(savename+'.png')

    return sp

def splat_tspec(filename,specnum=0,**kwargs):
    """
    Same as splat_1d for tspec data
    """

    tdata = pyfits.getdata(filename)
    theader = pyfits.getheader(filename)
    if len(tdata.shape) == 3:
        tdata = tdata[specnum,:,:]
    wavelength = tdata[0,:]
    spectrum   = tdata[1,:]
    error      = tdata[2,:]

    vconv = lambda x: wavelength[x]
    ivconv = lambda x: argmin(abs(wavelength-x))
    
    specname='TSPEC'
    dv = median(wavelength[1:] - wavelength[:-1])
    
    sp = SpecPlotter(spectrum,vconv=vconv,specname=specname,dv=dv,hdr=theader)

    sp.plotspec(cube=False,units=theader.get('YUNITS'),xunits=theader.get('XUNITS'),**kwargs)

    return sp
