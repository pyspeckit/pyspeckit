import matplotlib
import matplotlib.pyplot
import matplotlib.figure

class Plotter(object):
    """
    Class to plot a spectrum
    """


    def __init__(self,Spectrum,autorefresh=True):
        self.figure = None
        self.axis = None
        self.Spectrum = Spectrum

        # plot parameters
        self.offset = 0.0 # vertical offset
        self.autorefresh = autorefresh
        self.xlabel = ""
        self.ylabel = ""
        self.title  = ""

    def __call__(self, figure=None, axis=None, clear=True , **kwargs):
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
        elif isinstance(axis,matplotlib.axes.Axes):
            self.axis = axis
            self.figure = axis.figure
        elif type(figure) is int:
            self.figure = matplotlib.pyplot.figure(figure)
        else:
            self.figure = matplotlib.pyplot.figure()
        if len(self.figure.axes) > 0 and self.axis is None:
            self.axis = self.figure.axes[0] # default to first axis
        else:
            self.axis = self.figure.gca()

        if clear: self.axis.clear()

        self.plot(**kwargs)

    def plot(self, offset=0.0, color='k', linestyle='steps-mid', linewidth=0.5,
            vmin=None, vmax=None, **kwargs):

        self.offset += offset

        self._spectrumplot = self.axis.plot(self.Spectrum.xarr,
                self.Spectrum.data+self.offset, color=color,
                linestyle=linestyle, linewidth=linewidth, **kwargs)

        if vmin is not None: xlo = self.vmin
        else: xlo=self.Spectrum.xarr.min()
        if vmax is not None: xhi = self.vmax
        else: xhi=self.Spectrum.xarr.max()
        self.axis.set_xlim(xlo,xhi)
        
        self.label()
        if self.autorefresh: self.refresh()

    def label(self, title=None, xlabel=None):
   
        if title is not None:
            self.title = title
        if self.title is not "":
            self.axis.set_title(self.title)

        if xlabel is not None:
            self.xlabel = xlabel
        elif isinstance(self.Spectrum.xtype,str):
            self.xlabel = self.Spectrum.xtype
            if isinstance(self.Spectrum.xunits,str):
                self.xlabel += " "+self.Spectrum.xunits
        if self.xlabel is not None:
            self.axis.set_xlabel(self.xlabel)

        if self.Spectrum.units in ['Ta*','Tastar','K']:
          self.axis.set_ylabel("$T_A^*$ (K)")
        elif self.Spectrum.units == 'mJy':
          self.axis.set_ylabel("$S_\\nu$ (mJy)")
        elif self.Spectrum.units == 'Jy':
          self.axis.set_ylabel("$S_\\nu$ (Jy)")
        else:
          self.axis.set_ylabel(self.Spectrum.units)

    def refresh(self):
        self.axis.figure.canvas.draw()
