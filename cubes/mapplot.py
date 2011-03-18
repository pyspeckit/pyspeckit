import matplotlib
import matplotlib.pyplot
import matplotlib.figure
import numpy as np
import pyfits
import cubes
try:
    import aplpy
    icanhasaplpy = True
except ImportError:
    icanhasaplpy = False

class MapPlotter(object):
    """
    Class to plot a spectrum
    """

    def __init__(self, Cube=None, figure=None, **kwargs):
        """
        Create a map figure for future plotting
        """
        # figure out where to put the plot
        if isinstance(figure,matplotlib.figure.Figure):
            self.figure = figure
        elif type(figure) is int:
            self.figure = matplotlib.pyplot.figure(figure)
        else:
            self.figure = matplotlib.pyplot.figure()
        self.axis = None

        self.Cube = Cube

        self.mapplot(**kwargs)

    def mapplot(self, estimator=np.mean, **kwargs):
        """
        Plot up a map based on an input data cube
        """
        # THIS IS A HACK!!!  isinstance(a function, function) must be a thing...
        FUNCTION = type(np.max)
        if type(estimator) is FUNCTION:
            plane = estimator(self.Cube.cube,axis=0)
        elif type(estimator) is str:
            if estimator == 'max':
                plane = self.Cube.cube.max(axis=0)
            elif estimator == 'int':
                dx = np.abs(self.cube.xarr[1:] - self.cube.xarr[:-1])
                dx = np.concatenate(dx,dx[-1])
                plane = self.Cube.cube.sum(axis=0) * dx
            elif estimator[-5:] == ".fits":
                plane = pyfits.getdata(estimator)

        if icanhasaplpy:
            self.figure.clf()
            self.fitsfile = pyfits.PrimaryHDU(data=plane,header=cubes.flatten_header(self.Cube.header))
            self.FITSFigure = aplpy.FITSFigure(self.fitsfile,figure=self.figure)
            self.FITSFigure.show_colorscale()
            self.axis = self.FITSFigure._ax1
        else:
            if self.axis is None:
                self.axis = self.figure.add_subplot(111)
            self.axis.imshow(plane)

        self.canvas = self.axis.figure.canvas

        self.clickid = self.canvas.mpl_connect('button_press_event',self.plot_spectrum)

    def plot_spectrum(self, event):
        """
        Connects map cube to Spectrum...
        """
        if event.inaxes:
          clickX = event.xdata
          clickY = event.ydata

          # grab toolbar info so that we don't do anything if a tool is selected
          tb = self.canvas.toolbar
          if tb.mode != '':
              return
          elif event.button==1:
              print "Plotting spectrum from point %i,%i" % (clickX,clickY)
              self.Cube.plot_spectrum(clickY,clickX,clear=True)
          elif event.button==2:
              print "OverPlotting spectrum from point %i,%i" % (clickX,clickY)
              self.Cube.plot_spectrum(clickY,clickX,clear=False)
          elif event.button==3:
              print "Disconnecting GAIA-like tool"
              self.canvas.mpl_disconnect(self.clickid)
          else:
              print "Call failed for some reason: "
              print "event: ",event
