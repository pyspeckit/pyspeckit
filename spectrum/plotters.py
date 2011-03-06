import matplotlib

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
            if len(self.figure.axes) > 0:
                self.axis = self.figure.axes
        elif isinstance(axis,matplotlib.axes.Axes):
            self.axis = axis
            self.figure = axis.figure
        elif type(figure) is int:
            self.figure = matplotlib.pyplot.figure(figure)
        else:
            self.figure = matplotlib.pyplot.figure()

        if clear: self.axis.clear()

        self.plot(**kwargs)

    def plot(self, offset=0.0, color='k', linestyle='steps-mid', linewidth=0.5, **kwargs):

        self.offset += offset

        self._spectrumplot = self.axis.plot(self.Spectrum.xarr,
                self.Spectrum.data+self.offset, color=color,
                linestyle=linestyle, linewidth=linewidth, **kwargs)

        if self.autorefresh: self.refresh()

    def refresh(self):
        self.axis.figure.canvas.draw()
