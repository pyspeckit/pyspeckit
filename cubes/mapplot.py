"""
Author: Adam Ginsburg
Created: 03/17/2011
"""
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

    def __init__(self, Cube=None, figure=None, doplot=False, **kwargs):
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
        self.FITSFigure = None
        self._click_marks = []
        self._circles = []

        self.Cube = Cube

        if doplot: self.mapplot(**kwargs)

    def __call__(self, **kwargs):
        return self.mapplot(**kwargs)

    def mapplot(self, estimator=np.mean, convention='calabretta', **kwargs):
        """
        Plot up a map based on an input data cube
        """
        # THIS IS A HACK!!!  isinstance(a function, function) must be a thing...
        FUNCTION = type(np.max)
        if type(estimator) is FUNCTION:
            self.plane = estimator(self.Cube.cube,axis=0)
        elif type(estimator) is str:
            if estimator == 'max':
                self.plane = self.Cube.cube.max(axis=0)
            elif estimator == 'int':
                dx = np.abs(self.cube.xarr[1:] - self.cube.xarr[:-1])
                dx = np.concatenate(dx,dx[-1])
                self.plane = self.Cube.cube.sum(axis=0) * dx
            elif estimator[-5:] == ".fits":
                self.plane = pyfits.getdata(estimator)

        if icanhasaplpy:
            self.figure.clf()
            self.fitsfile = pyfits.PrimaryHDU(data=self.plane,header=cubes.flatten_header(self.Cube.header))
            vmin = self.plane[self.plane==self.plane].min()
            vmax = self.plane[self.plane==self.plane].max()
            self.FITSFigure = aplpy.FITSFigure(self.fitsfile,figure=self.figure,convention=convention)
            self.FITSFigure.show_colorscale(vmin=vmin,vmax=vmax)
            self.axis = self.FITSFigure._ax1
        else:
            if self.axis is None:
                self.axis = self.figure.add_subplot(111)
            self.axis.imshow(self.plane)

        self.canvas = self.axis.figure.canvas

        self.clickid = self.canvas.mpl_connect('button_press_event',self.plot_spectrum)
        self.keyid = self.canvas.mpl_connect('key_press_event',self.plot_spectrum)

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
          elif event.key is not None: #hasattr(event,'key'):
              if event.key == 'c':
                  self.center = (clickX,clickY)
                  self._remove_circle()
                  self._add_click_mark(clickX,clickY,clear=True)
              elif event.key == 'r':
                  x,y = self.center
                  self._add_circle(x,y,clickX,clickY)
                  self.circle(x,y,clickX,clickY)
          elif event.button is not None: #hasattr(event,'button'):
              if event.button==1:
                  print "Plotting spectrum from point %i,%i" % (clickX,clickY)
                  self._remove_circle()
                  self._add_click_mark(clickX,clickY,clear=True)
                  self.Cube.plot_spectrum(clickY,clickX,clear=True)
              elif event.button==2:
                  print "OverPlotting spectrum from point %i,%i" % (clickX,clickY)
                  self._add_click_mark(clickX,clickY,clear=False)
                  self.Cube.plot_spectrum(clickY,clickX,clear=False)
              elif event.button==3:
                  print "Disconnecting GAIA-like tool"
                  self.canvas.mpl_disconnect(self.clickid)
                  self.canvas.mpl_disconnect(self.keyid)
          else:
              print "Call failed for some reason: "
              print "event: ",event

    def _add_click_mark(self,x,y,clear=False):
        """
        Add an X at some position
        """
        if clear:
            self._clear_click_marks()
        if self.FITSFigure is not None:
            label = 'xmark%i' % (len(self._click_marks)+1)
            x,y = self.FITSFigure.pixel2world(x,y)
            self.FITSFigure.show_markers(x,y,marker='x',c='k',layer=label)
            self._click_marks.append( label )
        else:
            self._click_marks.append( self.axis.plot(x,y,'kx') )
        self.refresh()

    def _clear_click_marks(self):
        """
        Remove all marks added by previous clicks
        """
        if self.FITSFigure is not None:
            for mark in self._click_marks:
                if mark in self.FITSFigure._layers:
                    self.FITSFigure.remove_layer(mark)
        else:
            for mark in self._click_marks:
                self._click_marks.remove(mark)
                self.axis.lines.remove(mark)
            self.refresh()

    def _add_circle(self,x,y,x2,y2):
        """
        """
        if self.FITSFigure is not None:
            x,y = self.FITSFigure.pixel2world(x,y)
            x2,y2 = self.FITSFigure.pixel2world(x2,y2)
            r = (3600.0*np.linalg.norm(np.array([x,y])-np.array([x2,y2])))**2
            self.FITSFigure.show_markers(x,y,s=r,marker='o',facecolor='none',edgecolor='black',layer='circle')
        else:
            r = np.linalg.norm(np.array([x,y])-np.array([x2,y2]))
            self._circles.append( matplotlib.patches.Circle([x,y],radius=r) )
            self.axis.patches.append(circle)
            self.refresh()

    def _remove_circle(self):
        """
        """
        if self.FITSFigure is not None:
            if 'circle' in self.FITSFigure._layers:
                self.FITSFigure.remove_layer('circle')
        else:
            for circle in self._circles:
                self.axis.patches.remove(circle)
            self.refresh()

    def refresh(self):
        if self.axis is not None:
            self.axis.figure.canvas.draw()

    def circle(self,x1,y1,x2,y2):
        """
        Plot the spectrum of a circular aperture
        """

        r = (np.linalg.norm(np.array([x1,y1])-np.array([x2,y2])))
        self.Cube.plot_apspec([x1,y1,r])
        #self.Cube.data = cubes.extract_aperture( self.Cube.cube, [x1,y1,r] , coordsys=None )
        #self.Cube.plotter()
