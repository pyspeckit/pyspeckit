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

def Property(func):
    return property(**func())

def add_plot_property(propname,doc=""):

    def fget(self):
        return self.__dict__["_"+propname]

    def fset(self,new_prop):
        self.__dict__["_"+propname] = new_prop

    fdel = None

    return_dict = locals()
    del return_dict['propname']

    #import pdb; pdb.set_trace()
    #return property(get_prop,set_prop,None,prop_details)
    return return_dict

class Plotter(object):
    """
    Class to plot a spectrum
    """


    def __init__(self, Spectrum, autorefresh=True, title="", ylabel="",
            xlabel="", **kwargs):
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
        self._xmax = None
        self._xmin = None
        self._ymax = None
        self._ymin = None

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
        else:
            self.axis = self.figure.gca()


        if clear: self.axis.clear()

        self.plotkwargs = kwargs

        self.plot(**kwargs)

    @Property
    def ymax():
        return add_plot_property('ymax','Maximum Y value on plot')
        #self.ymax = add_plot_property(self,'ymax','Maximum Y value on plot')
        #self.ymin = add_plot_property(self, 
    @Property
    def ymin():
        return add_plot_property('ymin','Minimum Y value on plot')
        #self.xmax = add_plot_property(self,    
    @Property
    def xmax():
        return add_plot_property('xmax','Maximum X value on plot')
        #self.xmin = add_plot_property(self,    
    @Property
    def xmin():
        return add_plot_property('xmin','Minimum X value on plot')

    #@property
    #def get_ymax(self):
    #    return self._ymax

    #@ymax.setter
    #def set_ymax(self,new_ymax):
    #    self._ymax = new_ymax

    #@property
    #def get_ymax(self):
    #    return self._ymax

    #@ymax.setter
    #def set_ymax(self,new_ymax):
    #    self._ymax = new_ymax

    def plot(self, offset=0.0, color='k', linestyle='steps-mid', linewidth=0.5,
            errstyle=None, erralpha=0.2, silent=False, **kwargs):
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
        for arg in ['title','xlabel','ylabel']:
            if kwargs.has_key(arg): kwargs.pop(arg)

        reset_kwargs = {}
        for arg in ['xmin','xmax','ymin','ymax','reset_xlimits','reset_ylimits','ypeakscale']:
            if kwargs.has_key(arg): reset_kwargs[arg] = kwargs.pop(arg)

        self._spectrumplot = self.axis.plot(self.Spectrum.xarr,
                self.Spectrum.data+self.offset, color=color,
                linestyle=linestyle, linewidth=linewidth, **kwargs)
        if errstyle is not None:
            if errstyle == 'fill':
                self.errorplot = [self.axis.fill_between(steppify(self.Spectrum.xarr,isX=True),
                        steppify(self.Spectrum.data+self.offset-self.Spectrum.error),
                        steppify(self.Spectrum.data+self.offset+self.Spectrum.error),
                        facecolor=color, alpha=erralpha, **kwargs)]
            elif errstyle == 'bars':
                self.errorplot = axis.errorbar(self.Spectrum.xarr, self.Spectrum.data+self.offset,
                        yerr=self.Spectrum.error, ecolor=color, fmt=None,
                        **kwargs)

        self.reset_limits(silent=silent, **reset_kwargs)

        if self.autorefresh: self.refresh()
    
    def reset_limits(self,xmin=None, xmax=None, ymin=None, ymax=None,
            reset_xlimits=False, reset_ylimits=False, ypeakscale=1.2,
            silent=False, **kwargs):
        """
        Automatically or manually reset the plot limits
        """

        if (self.Spectrum.xarr.max() < self._xmin or self.Spectrum.xarr.min() > self._xmax 
                or reset_xlimits):
            if not silent: print "Resetting X-axis min/max because the plot is out of bounds."
            self._xmin = None
            self._xmax = None
        if xmin is not None: self._xmin = xmin
        elif self._xmin is None: self._xmin=self.Spectrum.xarr.min()
        if xmax is not None: self._xmax = xmax
        elif self._xmax is None: self._xmax=self.Spectrum.xarr.max()
        self.axis.set_xlim(self._xmin,self._xmax)

        xpixmin = np.argmin(np.abs(self.Spectrum.xarr-self._xmin))
        xpixmax = np.argmin(np.abs(self.Spectrum.xarr-self._xmax))
        if xpixmin>xpixmax: xpixmin,xpixmax = xpixmax,xpixmin
        
        if (self.Spectrum.data.max() < self._ymin or self.Spectrum.data.min() > self._ymax
                or reset_ylimits):
            if not silent: print "Resetting Y-axis min/max because the plot is out of bounds."
            self._ymin = None
            self._ymax = None
        if ymin is not None: self._ymin = ymin
        elif self._ymin is None: self._ymin=np.nanmin(self.Spectrum.data[xpixmin:xpixmax])
        if ymax is not None: self._ymax = ymax
        elif self._ymax is None: self._ymax=(np.nanmax(self.Spectrum.data[xpixmin:xpixmax])-self._ymin) * ypeakscale + self._ymin
        self.axis.set_ylim(self._ymin,self._ymax)
        

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
