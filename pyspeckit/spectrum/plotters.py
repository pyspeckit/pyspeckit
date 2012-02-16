"""
=======
Plotter
=======

.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>
"""
import matplotlib
import matplotlib.pyplot
import matplotlib.figure
import itertools
from ..config import *
import numpy as np
from pyspeckit.specwarnings import warn
import copy

interactive_help_message = """
Interactive key commands for plotter.  An additional help message may appear if
you have initiated the fitter.
'?' - bring up this message
'f' - initiate the fitter 
'b' - initiate the baseliner 
'B' - initiate the baseliner (reset the selection too)
"""

class Plotter(object):
    """
    Class to plot a spectrum
    """


    def __init__(self, Spectrum, autorefresh=True, title="", ylabel="",
            xlabel="", silent=False, plotscale=1.0, **kwargs):
        self.figure = None
        self.axis = None
        self.Spectrum = Spectrum

        # plot parameters
        self.offset = 0.0 # vertical offset
        self.autorefresh = autorefresh
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.title  = title
        self.errorplot = None
        self.plotkwargs = kwargs
        self.xmax = None
        self.xmin = None
        self.ymax = None
        self.ymin = None
        self.keyclick = None
        self.silent = silent
        self.plotscale = plotscale

        self._xclick1 = None
        self._xclick2 = None

    def __call__(self, figure=None, axis=None, clear=True, autorefresh=None, plotscale=1.0, **kwargs):
        """
        Plot a spectrum
        
        Keywords:
        figure - either a matplotlib figure instance or a figure number
            to pass into pyplot.figure.  
        axis - Alternative to figure, can pass an axis instance and use
            it as the plotting canvas
        clear - Clear the axis before plotting?
        """

        # figure out where to put the plot
        if isinstance(figure,matplotlib.figure.Figure):
            self.figure = figure
            self.axis = self.figure.gca()
        elif type(figure) is int:
            self.figure = matplotlib.pyplot.figure(figure)
            self.axis = self.figure.gca()
        elif self.figure is None:
            if isinstance(axis,matplotlib.axes.Axes):
                self.axis = axis
                self.figure = axis.figure
            else:
                self.figure = matplotlib.pyplot.figure()

        if self.keyclick is None:
            self.keyclick = self.figure.canvas.mpl_connect('key_press_event',self.parse_keys)

        if axis is not None:
            self.figure.canvas.mpl_disconnect(self.keyclick)
            self.axis = axis
            self.figure = axis.figure
            self.keyclick = self.figure.canvas.mpl_connect('key_press_event',self.parse_keys)
        elif len(self.figure.axes) > 0 and self.axis is None:
            self.axis = self.figure.axes[0] # default to first axis
        elif self.axis is None:
            self.axis = self.figure.gca()


        if clear and self.axis is not None: self.axis.clear()

        if autorefresh is not None:
            self.autorefresh = autorefresh

        self.plotscale = plotscale

        self.plotkwargs = kwargs

        self.plot(**kwargs)

    def plot(self, offset=0.0, xoffset=0.0, color='k', linestyle='steps-mid',
            linewidth=0.5, errstyle=None, erralpha=0.2, silent=None,
            reset=True, refresh=True, **kwargs):
        """
        Plot the spectrum!

        Tries to automatically find a reasonable plotting range if one is not set.  

        offset - vertical offset to add to the spectrum before plotting.  Useful if you
        want to overlay multiple spectra on a single plot

        color - default to plotting spectrum in black

        linestyle - histogram-style plotting
        linewidth - narrow lines are helpful when histo-plotting
        xmin/xmax/ymin/ymax - override defaults for plot range.  Once set, these parameters
        are sticky (i.e., replotting will use the same ranges)

        reset_[xy]limits - Reset the limits to "sensible defaults" 

        ypeakscale - Scale up the Y maximum value.  Useful to keep the
        annotations away from the data.

        errstyle - can be "fill", which draws partially transparent boxes around the data to show
            the error region, or "bars" which draws standard errorbars

        *reset* [ bool ]
            Reset the x/y axis limits?
        """

        if self.axis is None:
            raise Exception("You must call the Plotter class to initiate the canvas before plotting.")

        self.offset += offset

        self.label(**kwargs)
        for arg in ['title','xlabel','ylabel']:
            if arg in kwargs: kwargs.pop(arg)

        reset_kwargs = {}
        for arg in ['xmin','xmax','ymin','ymax','reset_xlimits','reset_ylimits','ypeakscale']:
            if arg in kwargs: reset_kwargs[arg] = kwargs.pop(arg)

        self._spectrumplot = self.axis.plot(self.Spectrum.xarr+xoffset,
                self.Spectrum.data*self.plotscale+self.offset, color=color,
                linestyle=linestyle, linewidth=linewidth, **kwargs)

        if errstyle is not None:
            if errstyle == 'fill':
                order = -1 if self.Spectrum.xarr[-1] < self.Spectrum.xarr[0] else 1
                self.errorplot = [self.axis.fill_between(steppify(self.Spectrum.xarr[::order]+xoffset,isX=True),
                    steppify((self.Spectrum.data*self.plotscale+self.offset-self.Spectrum.error*self.plotscale)[::order]),
                    steppify((self.Spectrum.data*self.plotscale+self.offset+self.Spectrum.error*self.plotscale)[::order]),
                    facecolor=color, alpha=erralpha, **kwargs)]
            elif errstyle == 'bars':
                self.errorplot = self.axis.errorbar(self.Spectrum.xarr+xoffset, self.Spectrum.data*self.plotscale+self.offset,
                        yerr=self.Spectrum.error*self.plotscale, ecolor=color, fmt=None,
                        **kwargs)

        if silent is not None:
            self.silent = silent

        if reset:
            self.reset_limits(**reset_kwargs)

        if self.autorefresh and refresh: self.refresh()
    
    def reset_limits(self,xmin=None, xmax=None, ymin=None, ymax=None,
            reset_xlimits=False, reset_ylimits=False, ypeakscale=1.2,
            silent=None, **kwargs):
        """
        Automatically or manually reset the plot limits
        """

        if silent is not None:
            self.silent = silent

        if (self.Spectrum.xarr.max() < self.xmin or self.Spectrum.xarr.min() > self.xmax 
                or reset_xlimits):
            if not self.silent: warn( "Resetting X-axis min/max because the plot is out of bounds." )
            self.xmin = None
            self.xmax = None
        if xmin is not None: self.xmin = xmin
        elif self.xmin is None: self.xmin=self.Spectrum.xarr.min()
        if xmax is not None: self.xmax = xmax
        elif self.xmax is None: self.xmax=self.Spectrum.xarr.max()
        self.axis.set_xlim(self.xmin,self.xmax)

        xpixmin = np.argmin(np.abs(self.Spectrum.xarr-self.xmin))
        xpixmax = np.argmin(np.abs(self.Spectrum.xarr-self.xmax))
        if xpixmin>xpixmax: xpixmin,xpixmax = xpixmax,xpixmin
        elif xpixmin == xpixmax:
            if not self.silent: warn( "ERROR: the X axis limits specified were invalid.  Resetting." )
            self.reset_limits(reset_xlimits=True, ymin=ymin, ymax=ymax, reset_ylimits=reset_ylimits, ypeakscale=ypeakscale, **kwargs)
            return
        
        if (self.Spectrum.data.max() < self.ymin or self.Spectrum.data.min() > self.ymax
                or reset_ylimits):
            if not self.silent and not reset_ylimits: warn( "Resetting Y-axis min/max because the plot is out of bounds." )
            self.ymin = None
            self.ymax = None
        if ymin is not None: self.ymin = ymin
        elif self.ymin is None: self.ymin=np.nanmin(self.Spectrum.data[xpixmin:xpixmax],axis=0)
        if ymax is not None: self.ymax = ymax
        elif self.ymax is None: self.ymax=(np.nanmax(self.Spectrum.data[xpixmin:xpixmax],axis=0)-self.ymin) * ypeakscale + self.ymin
        self.axis.set_ylim(self.ymin+self.offset,self.ymax+self.offset)
        

    def label(self, title=None, xlabel=None, ylabel=None, **kwargs):
   
        if title is not None:
            self.title = title
        elif hasattr(self.Spectrum,'specname'):
            self.title = self.Spectrum.specname
        if self.title is not "":
            self.axis.set_title(self.title)

        if xlabel is not None:
            self.xlabel = xlabel
        elif isinstance(self.Spectrum.xarr.xtype,str):
            self.xlabel = self.Spectrum.xarr.xtype.title()
            if isinstance(self.Spectrum.xarr.units,str):
                self.xlabel += " ("+self.Spectrum.xarr.units+")"
        if self.xlabel is not None:
            self.axis.set_xlabel(self.xlabel)

        if ylabel is not None:
            self.ylabel=ylabel
            self.axis.set_ylabel(self.ylabel)
        elif self.Spectrum.units in ['Ta*','Tastar','K']:
            self.axis.set_ylabel("$T_A^*$ (K)")
        elif self.Spectrum.units == 'mJy':
            self.axis.set_ylabel("$S_\\nu$ (mJy)")
        elif self.Spectrum.units == 'Jy':
            self.axis.set_ylabel("$S_\\nu$ (Jy)")            
        else:
            if "$" in self.Spectrum.units:
                # assume LaTeX already
                self.axis.set_ylabel(self.Spectrum.units)
            elif len(self.Spectrum.units.split()) > 1: 
                if self.Spectrum.units.split()[1] in ['erg/cm^2/s/Ang',
                        'erg/cm2/s/A', 'erg/cm2/s/Ang', 'erg/cm/s/Ang']:
                    norm = parse_norm(self.Spectrum.units.split()[0])
                    self.axis.set_ylabel("$%s \\mathrm{erg/s/cm^2/\\AA}$" % norm)
                elif self.Spectrum.units.split()[1] in ['W/m^2/Hz','w/m^2/hz','W/m/hz','W/m/Hz']:
                    norm = parse_norm(self.Spectrum.units.split()[0])
                    self.axis.set_ylabel("$%s \\mathrm{W/m^2/Hz}$" % norm)
            else:
                label_units = parse_units(self.Spectrum.units)
                self.axis.set_ylabel(label_units)

    def refresh(self):
        if self.axis is not None:
            self.axis.figure.canvas.draw()

    def savefig(self,fname,bbox_inches='tight',**kwargs):
        """
        simple wrapper of maplotlib's savefig.  
        """
        self.axis.figure.savefig(fname,bbox_inches=bbox_inches,**kwargs)

    def parse_keys(self,event):
        """
        Parse key commands entered from the keyboard
        """
        if hasattr(event,'key'):
            if event.key == '?':
                print interactive_help_message
            elif event.key == 'f':
                print "Fitter initiated from the interactive plotter."
                self.Spectrum.specfit(interactive=True)
            elif event.key == 'b':
                print "Baseline initiated from the interactive plotter."
                self.Spectrum.baseline(interactive=True, reset_selection=False)
            elif event.key == 'B':
                print "Baseline initiated from the interactive plotter (with reset)."
                self.Spectrum.baseline(interactive=True, reset_selection=True)

    def get_two_clicks(self,event):

        if self._xclick1 is None:
            self._xclick1 = event.xdata
        elif self._xclick2 is None:
            self._xclick2 = event.xdata

    def copy(self, parent=None):
        """
        Create a copy of the plotter with blank (uninitialized) axis & figure

        [ parent ] 
            A spectroscopic axis instance that is the parent of the specfit
            instance.  This needs to be specified at some point, but defaults
            to None to prevent overwriting a previous plot.
        """

        newplotter = copy.copy(self)
        newplotter.Spectrum = parent
        newplotter.axis = None
        newplotter.figure = None

        return newplotter

def parse_units(labelstring):
    import re
    labelstring = re.sub("um","$\mu$m",labelstring)
    labelstring = re.sub("-1","$^{-1}$",labelstring)
    labelstring = re.sub("-2","$^{-2}$",labelstring)
    labelstring = re.sub("-3","$^{-3}$",labelstring)
    labelstring = re.sub("ergss","ergs s",labelstring)
    return labelstring
    
def parse_norm(norm):
    """
    Expected format: norm = 10E15
    """    
        
    try: 
        base, exp = norm.split('E')
    except ValueError:
        base, exp = norm.split('e')

    if float(base) == 1.0:
        norm = '10'
    else:
        norm = base    
        
    norm += '^{%s}' % exp    
    
    return norm

def steppify(arr,isX=False):
    """
    *support function*
    Converts an array to double-length for step plotting
    """
    if isX:
        interval = abs(arr[1:]-arr[:-1]) / 2.0
        newarr = np.array(zip(arr[:-1]-interval,arr[:-1]+interval)).ravel()
        newarr = np.concatenate([newarr,2*[newarr[-1]+interval[-1]]])
    else:
        newarr = np.array(zip(arr,arr)).ravel()
    return newarr
