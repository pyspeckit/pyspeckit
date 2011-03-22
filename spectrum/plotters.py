import matplotlib
import matplotlib.pyplot
import matplotlib.figure
from config import *
import numpy as np

interactive_help_message = """
Interactive key commands for plotter.  An additional help message may appear if
you have initiated the fitter.
'?' - bring up this message
'f' - initiate the fitter 
'b' - initiate the baseliner 
"""

class Plotter(object):
    """
    Class to plot a spectrum
    """


    def __init__(self, Spectrum, autorefresh=True, title="", ylabel="",
            xlabel=""):
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
        self.plotkwargs = {}
        self.xmax = None
        self.xmin = None
        self.ymax = None
        self.ymin = None
        self.keyclick = None

    def __call__(self, figure=None, axis=None, clear=True, **kwargs):
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
        if self.figure is None:
            if isinstance(figure,matplotlib.figure.Figure):
                self.figure = figure
            elif isinstance(axis,matplotlib.axes.Axes):
                self.axis = axis
                self.figure = axis.figure
            elif type(figure) is int:
                self.figure = matplotlib.pyplot.figure(figure)
            else:
                self.figure = matplotlib.pyplot.figure()

        if self.keyclick is None:
            self.keyclick = self.figure.canvas.mpl_connect('key_press_event',self.parse_keys)

        if len(self.figure.axes) > 0 and self.axis is None:
            self.axis = self.figure.axes[0] # default to first axis
        elif axis is not None:
            self.figure.canvas.disconnect(self.keyclick)
            self.axis = axis
            self.figure = axis.figure
            self.keyclick = self.figure.canvas.mpl_connect('key_press_event',self.parse_keys)
        else:
            self.axis = self.figure.gca()


        if clear: self.axis.clear()

        self.plot(**kwargs)

    def plot(self, offset=0.0, color='k', linestyle='steps-mid', linewidth=0.5,
            xmin=None, xmax=None, ymin=None, ymax=None, reset_xlimits=False,
            reset_ylimits=False, ypeakscale=1.2, **kwargs):
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
        """

        if self.axis is None:
            raise Exception("You must call the Plotter class to initiate the canvas before plotting.")

        self.offset += offset

        self.label(**kwargs)
        if kwargs.has_key('title'): kwargs.pop('title')
        if kwargs.has_key('xlabel'): kwargs.pop('xlabel')
        if kwargs.has_key('ylabel'): kwargs.pop('ylabel')
        self._spectrumplot = self.axis.plot(self.Spectrum.xarr,
                self.Spectrum.data+self.offset, color=color,
                linestyle=linestyle, linewidth=linewidth, **kwargs)

        if (self.Spectrum.xarr.max() < self.xmin or self.Spectrum.xarr.min() > self.xmax 
                or reset_xlimits):
            print "Resetting X-axis min/max because the plot is out of bounds."
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
        
        if (self.Spectrum.data.max() < self.ymin or self.Spectrum.data.min() > self.ymax
                or reset_ylimits):
            print "Resetting Y-axis min/max because the plot is out of bounds."
            self.ymin = None
            self.ymax = None
        if ymin is not None: self.ymin = ymin
        elif self.ymin is None: self.ymin=self.Spectrum.data[xpixmin:xpixmax].min()
        if ymax is not None: self.ymax = ymax
        elif self.ymax is None: self.ymax=self.Spectrum.data[xpixmin:xpixmax].max() * ypeakscale
        self.axis.set_ylim(self.ymin,self.ymax)
        
        if self.autorefresh: self.refresh()

    def label(self, title=None, xlabel=None, ylabel=None):
   
        if title is not None:
            self.title = title
        if self.title is not "":
            self.axis.set_title(self.title)

        if xlabel is not None:
            self.xlabel = xlabel
        elif isinstance(self.Spectrum.xarr.xtype,str):
            self.xlabel = self.Spectrum.xarr.xtype
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
                self.Spectrum.specfit(interactive=True)
            elif event.key == 'b':
                self.Spectrum.baseline(interactive=True)

def parse_units(labelstring):
    import re
    labelstring = re.sub("um","$\mu$m",labelstring)
    labelstring = re.sub("-1","$^{-1}$",labelstring)
    labelstring = re.sub("-2","$^{-2}$",labelstring)
    labelstring = re.sub("-3","$^{-3}$",labelstring)
    labelstring = re.sub("ergss","ergs s",labelstring)
    return labelstring
