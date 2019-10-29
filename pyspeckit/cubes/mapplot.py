"""
MapPlot
-------

Make plots of the cube and interactively connect them to spectrum plotting.
This is really an interactive component of the package; nothing in here is
meant for publication-quality plots, but more for user interactive analysis.

That said, the plotter makes use of `APLpy <https://github.com/aplpy/aplpy>`_,
so it is possible to make publication-quality plots.

:author: Adam Ginsburg
:date: 03/17/2011

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"""
from __future__ import print_function
import matplotlib
import matplotlib.figure
import numpy as np
import copy
import itertools
import six
try:
    import astropy.wcs as pywcs
    import astropy.io.fits as pyfits
    pywcsOK = True
except ImportError:
    try:
        import pyfits
        import pywcs
        pywcsOK = True
    except ImportError:
        pywcsOK = False
try:
    import aplpy
    icanhasaplpy = True
except: # aplpy fails with generic exceptions instead of ImportError
    icanhasaplpy = False

from . import cubes

class MapPlotter(object):
    """
    Class to plot a spectrum

    See `mapplot` for use documentation; this docstring is only for
    initialization.
    """

    def __init__(self, Cube=None, figure=None, doplot=False, **kwargs):
        """
        Create a map figure for future plotting
        """
        import matplotlib.pyplot
        self._pyplot = matplotlib.pyplot

        # figure out where to put the plot
        if isinstance(figure,matplotlib.figure.Figure):
            self.figure = figure
        elif type(figure) is int:
            self.figure = self._pyplot.figure(figure)
        else:
            self.figure = None
        self.axis = None
        self.FITSFigure = None
        self._click_marks = []
        self._circles = []
        self._clickX = None
        self._clickY = None

        self.overplot_colorcycle = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y'])
        self.overplot_linestyle = '-'

        self.Cube = Cube
        if self.Cube is not None:
            self.header = cubes.flatten_header(self.Cube.header, delete=True)
        if pywcsOK:
            self.wcs = pywcs.WCS(self.header)

        if doplot: self.mapplot(**kwargs)

    def __call__(self, **kwargs):
        """ see mapplot """
        return self.mapplot(**kwargs)

    def mapplot(self, convention='calabretta', colorbar=True, useaplpy=True,
                vmin=None, vmax=None, cmap=None, plotkwargs={}, **kwargs):
        """
        Plot up a map based on an input data cube.

        The map to be plotted is selected using `makeplane`.
        The `estimator` keyword argument is passed to that function.

        The plotted map, once shown, is interactive.  You can click on it with any
        of the three mouse buttons.

        Button 1 or keyboard '1':
            Plot the selected pixel's spectrum in another window.  Mark the
            clicked pixel with an 'x'
        Button 2 or keyboard 'o':
            Overplot a second (or third, fourth, fifth...) spectrum in the
            external plot window
        Button 3:
            Disconnect the interactive viewer

        You can also click-and-drag with button 1 to average over a circular
        region.  This same effect can be achieved by using the 'c' key to
        set the /c/enter of a circle and the 'r' key to set its /r/adius (i.e.,
        hover over the center and press 'c', then hover some distance away and
        press 'r').


        Parameters
        ----------
        convention : 'calabretta' or 'griesen'
            The default projection to assume for Galactic data when plotting
            with aplpy.
        colorbar : bool
            Whether to show a colorbar
        plotkwargs : dict, optional
            A dictionary of keyword arguments to pass to aplpy.show_colorscale
            or matplotlib.pyplot.imshow
        useaplpy : bool
            Use aplpy if a FITS header is available
        vmin, vmax: float or None
            Override values for the vmin/vmax values.  Will be automatically
            determined if left as None

        .. todo:
            Allow mapplot in subfigure
        """
        if (self.figure is None):
            self.figure = self._pyplot.figure()
        elif (not self._pyplot.fignum_exists(self.figure.number)):
            self.figure = self._pyplot.figure()
        else:
            self._disconnect()
            self.figure.clf()

        # this is where the map is created; everything below this is just plotting
        self.makeplane(**kwargs)

        # have tot pop out estimator so that kwargs can be passed to imshow
        if 'estimator' in kwargs:
            kwargs.pop('estimator')

        # Below here is all plotting stuff

        if vmin is None: vmin = self.plane[self.plane==self.plane].min()
        if vmax is None: vmax = self.plane[self.plane==self.plane].max()

        if icanhasaplpy and useaplpy:
            self.fitsfile = pyfits.PrimaryHDU(data=self.plane,header=self.header)
            self.FITSFigure = aplpy.FITSFigure(self.fitsfile,figure=self.figure,convention=convention)
            self.FITSFigure.show_colorscale(vmin=vmin, vmax=vmax, cmap=cmap, **plotkwargs)
            if hasattr(self.FITSFigure, '_ax1'):
                self.axis = self.FITSFigure._ax1
            else:
                self.axis = self.FITSFigure.ax
            if colorbar:
                try:
                    self.FITSFigure.add_colorbar()
                except Exception as ex:
                    print("ERROR: Could not create colorbar!  Error was %s" % str(ex))
            self._origin = 0 # FITS convention
            # TODO: set _origin to 1 if using PIXEL units, not real wcs
        else:
            self.axis = self.figure.add_subplot(111)
            if hasattr(self,'colorbar') and self.colorbar is not None:
                if self.colorbar.ax in self.axis.figure.axes:
                    self.axis.figure.delaxes(self.colorbar.ax)
            self.axis.imshow(self.plane, vmin=vmin, vmax=vmax, cmap=cmap, **plotkwargs)
            if colorbar:
                try:
                    self.colorbar = self._pyplot.colorbar(self.axis.images[0])
                except Exception as ex:
                    print("ERROR: Could not create colorbar!  Error was %s" % str(ex))
            self._origin = 0 # normal convention

        self.canvas = self.axis.figure.canvas

        self._connect()

    def _connect(self):
        """ Connect click, click up (release click), and key press to events """
        self.clickid = self.canvas.callbacks.connect('button_press_event',self.click)
        self.clickupid = self.canvas.callbacks.connect('button_release_event',self.plot_spectrum)
        self.keyid = self.canvas.callbacks.connect('key_press_event',self.plot_spectrum)

    def _disconnect(self):
        """ Disconnect click, click up (release click), and key press from events """
        if hasattr(self,'canvas'):
            self.canvas.mpl_disconnect(self.clickid)
            self.canvas.mpl_disconnect(self.clickupid)
            self.canvas.mpl_disconnect(self.keyid)

    def makeplane(self, estimator=np.nanmean):
        """
        Create a "plane" view of the cube, either by slicing or projecting it
        or by showing a slice from the best-fit model parameter cube.

        Parameters
        ----------

        estimator : [ function | 'max' | 'int' | FITS filename | integer | slice ]
            A non-pythonic, non-duck-typed variable.  If it's a function, apply that function
            along the cube's spectral axis to obtain an estimate (e.g., mean, min, max, etc.).
            'max' will do the same thing as passing np.max
            'int' will attempt to integrate the image (which is why I didn't duck-type)
            (integrate means sum and multiply by dx)
            a .fits filename will be read using pyfits (so you can make your own cover figure)
            an integer will get the n'th slice in the parcube if it exists
            If it's a slice, slice the input data cube along the Z-axis with this slice

        """
        # THIS IS A HACK!!!  isinstance(a function, function) must be a thing...
        FUNCTION = type(np.max)

        # estimator is NOT duck-typed
        if type(estimator) is FUNCTION:
            self.plane = estimator(self.Cube.cube,axis=0)
        elif isinstance(estimator, six.string_types):
            if estimator == 'max':
                self.plane = self.Cube.cube.max(axis=0)
            elif estimator == 'int':
                dx = np.abs(self.Cube.xarr[1:] - self.Cube.xarr[:-1])
                dx = np.concatenate([dx,[dx[-1]]])
                self.plane = (self.Cube.cube * dx[:,np.newaxis,np.newaxis]).sum(axis=0)
            elif estimator[-5:] == ".fits":
                self.plane = pyfits.getdata(estimator)
        elif type(estimator) is slice:
            self.plane = self.Cube.cube[estimator,:,:]
        elif type(estimator) is int:
            if hasattr(self.Cube,'parcube'):
                self.plane = self.Cube.parcube[estimator,:,:]

        if self.plane is None:
            raise ValueError("Invalid estimator %s" % (str(estimator)))

        if np.sum(np.isfinite(self.plane)) == 0:
            raise ValueError("Map is all NaNs or infs.  Check your estimator or your input cube.")

    def click(self,event):
        """
        Record location of downclick
        """
        if event.inaxes:
            self._clickX = np.round(event.xdata) - self._origin
            self._clickY = np.round(event.ydata) - self._origin

    def plot_spectrum(self, event, plot_fit=True):
        """
        Connects map cube to Spectrum...
        """
        self.event = event
        if event.inaxes:
            clickX = np.round(event.xdata) - self._origin
            clickY = np.round(event.ydata) - self._origin

            # grab toolbar info so that we don't do anything if a tool is selected
            tb = self.canvas.toolbar
            if tb.mode != '':
                return
            elif event.key is not None:
                if event.key == 'c':
                    self._center = (clickX-1,clickY-1)
                    self._remove_circle()
                    self._add_click_mark(clickX,clickY,clear=True)
                elif event.key == 'r':
                    x,y = self._center
                    self._add_circle(x,y,clickX,clickY)
                    self.circle(x,y,clickX-1,clickY-1)
                elif event.key == 'o':
                    clickX,clickY = round(clickX),round(clickY)
                    print("OverPlotting spectrum from point %i,%i" % (clickX-1,clickY-1))
                    color = next(self.overplot_colorcycle)
                    self._add_click_mark(clickX,clickY,clear=False, color=color)
                    self.Cube.plot_spectrum(clickX-1,clickY-1,clear=False, color=color, linestyle=self.overplot_linestyle)
                elif event.key in ('1','2'):
                    event.button = int(event.key)
                    event.key = None
                    self.plot_spectrum(event)
            elif (hasattr(event,'button') and event.button in (1,2)
                    and not (self._clickX == clickX and self._clickY == clickY)):
                if event.button == 1:
                    self._remove_circle()
                    clear=True
                    color = 'k'
                    linestyle = 'steps-mid'
                else:
                    color = next(self.overplot_colorcycle)
                    linestyle = self.overplot_linestyle
                    clear=False
                rad = ( (self._clickX-clickX)**2 + (self._clickY-clickY)**2 )**0.5
                print("Plotting circle from point %i,%i to %i,%i (r=%f)" % (self._clickX,self._clickY,clickX,clickY,rad))
                self._add_circle(self._clickX,self._clickY,clickX,clickY)
                self.circle(self._clickX,self._clickY,clickX,clickY,clear=clear,linestyle=linestyle,color=color)
            elif hasattr(event,'button') and event.button is not None:
                if event.button==1:
                    clickX,clickY = round(clickX),round(clickY)
                    print("Plotting spectrum from point %i,%i" % (clickX,clickY))
                    self._remove_circle()
                    self._add_click_mark(clickX,clickY,clear=True)
                    self.Cube.plot_spectrum(clickX,clickY,clear=True)
                    if plot_fit: self.Cube.plot_fit(clickX, clickY, silent=True)
                elif event.button==2:
                    clickX,clickY = round(clickX),round(clickY)
                    print("OverPlotting spectrum from point %i,%i" % (clickX,clickY))
                    color = next(self.overplot_colorcycle)
                    self._add_click_mark(clickX,clickY,clear=False, color=color)
                    self.Cube.plot_spectrum(clickX,clickY,clear=False, color=color, linestyle=self.overplot_linestyle)
                elif event.button==3:
                    print("Disconnecting GAIA-like tool")
                    self._disconnect()
            else:
                print("Call failed for some reason: ")
                print("event: ",event)
        else:
            pass
            # never really needed... warn("Click outside of axes")

    def _add_click_mark(self,x,y,clear=False,color='k'):
        """
        Add an X at some position
        """
        if clear:
            self._clear_click_marks()
        if self.FITSFigure is not None:
            label = 'xmark%i' % (len(self._click_marks)+1)
            x,y = self.FITSFigure.pixel2world(x,y)
            self.FITSFigure.show_markers(x,y,marker='x',c=color,layer=label)
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
                if mark in self.axis.lines:
                    self.axis.lines.remove(mark)
            self.refresh()

    def _add_circle(self,x,y,x2,y2,**kwargs):
        """
        """
        if self.FITSFigure is not None:
            x,y = self.FITSFigure.pixel2world(x,y)
            x2,y2 = self.FITSFigure.pixel2world(x2,y2)
            r = (np.linalg.norm(np.array([x,y])-np.array([x2,y2])))
            #self.FITSFigure.show_markers(x,y,s=r,marker='o',facecolor='none',edgecolor='black',layer='circle')
            layername = "circle%02i" % len(self._circles)
            self.FITSFigure.show_circles(x,y,r,edgecolor='black',facecolor='none',layer=layername,**kwargs)
            self._circles.append(layername)
        else:
            r = np.linalg.norm(np.array([x,y])-np.array([x2,y2]))
            circle = matplotlib.patches.Circle([x,y],radius=r,**kwargs)
            self._circles.append( circle )
            self.axis.patches.append(circle)
            self.refresh()

    def _remove_circle(self):
        """
        """
        if self.FITSFigure is not None:
            for layername in self._circles:
                if layername in self.FITSFigure._layers:
                    self.FITSFigure.remove_layer(layername)
        else:
            for circle in self._circles:
                if circle in self.axis.patches:
                    self.axis.patches.remove(circle)
                self._circles.remove(circle)
            self.refresh()

    def refresh(self):
        if self.axis is not None:
            self.axis.figure.canvas.draw()

    def circle(self,x1,y1,x2,y2,**kwargs):
        """
        Plot the spectrum of a circular aperture
        """

        r = (np.linalg.norm(np.array([x1,y1])-np.array([x2,y2])))
        self.Cube.plot_apspec([x1,y1,r],**kwargs)
        #self.Cube.data = cubes.extract_aperture( self.Cube.cube, [x1,y1,r] , coordsys=None )
        #self.Cube.plotter()

    def copy(self, parent=None):
        """
        Create a copy of the map plotter with blank (uninitialized) axis & figure

        [ parent ]
            A spectroscopic axis instance that is the parent of the specfit
            instance.  This needs to be specified at some point, but defaults
            to None to prevent overwriting a previous plot.
        """

        newmapplot = copy.copy(self)
        newmapplot.Cube = parent
        newmapplot.axis = None
        newmapplot.figure = None

        return newmapplot
