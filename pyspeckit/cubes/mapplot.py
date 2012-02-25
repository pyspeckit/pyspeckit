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
import matplotlib
import matplotlib.pyplot
import matplotlib.figure
import numpy as np
import pyfits
import copy
import itertools
from pyspeckit.specwarnings import warn
try:
    # pywcs is preferred over astropy.wcs because astropy.wcs failed on me
    import pywcs
    pywcsOK = True
except ImportError:
    try:
        import astropy.wcs as pywcs
        import astropy.io.fits as pyfits
        pywcsOK = True
    except ImportError:
        pywcsOK = False
import cubes
try:
    import aplpy
    icanhasaplpy = True
except: # aplpy fails with generic exceptions instead of ImportError 
    icanhasaplpy = False
try:
    import coords
    icanhascoords = True
except ImportError:
    icanhascoords = False


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
            self.header = cubes.flatten_header(self.Cube.header)
        if pywcsOK:
            self.wcs = pywcs.WCS(self.header)

        if doplot: self.mapplot(**kwargs)

    def __call__(self, **kwargs):
        """ see mapplot """
        return self.mapplot(**kwargs)

    def mapplot(self, convention='calabretta', colorbar=True, useaplpy=True,
            vmin=None, vmax=None, **kwargs):
        """
        Plot up a map based on an input data cube.

        The map to be plotted is selected using :function:`makeplane`.
        The :param:`estimator` keyword argument is passed to that function.

        `kwargs` are passed to aplpy.show_colorscale or
        matplotlib.pyplot.imshow (depending on whether aplpy is installed)

        TODO: Allow mapplot in subfigure
        """
        if self.figure is None:
            self.figure = matplotlib.pyplot.figure()

        self.makeplane(**kwargs)
        if 'estimator' in kwargs:
            kwargs.pop('estimator')

        if vmin is None: vmin = self.plane[self.plane==self.plane].min()
        if vmax is None: vmax = self.plane[self.plane==self.plane].max()

        if icanhasaplpy and useaplpy:
            self.figure.clf()
            self.fitsfile = pyfits.PrimaryHDU(data=self.plane,header=self.header)
            self.FITSFigure = aplpy.FITSFigure(self.fitsfile,figure=self.figure,convention=convention)
            self.FITSFigure.show_colorscale(vmin=vmin, vmax=vmax, **kwargs)
            self.axis = self.FITSFigure._ax1
            if colorbar: self.FITSFigure.add_colorbar()
        else:
            if self.axis is None:
                self.axis = self.figure.add_subplot(111)
            self.axis.imshow(self.plane, vmin=vmin, vmax=vmax, **kwargs)
            if colorbar: self.colorbar = matplotlib.pyplot.colorbar()

        self.canvas = self.axis.figure.canvas

        self._connect()

    def _connect(self):
        """ Connect click, click up (release click), and key press to events """
        self.clickid = self.canvas.mpl_connect('button_press_event',self.click)
        self.clickupid = self.canvas.mpl_connect('button_release_event',self.plot_spectrum)
        self.keyid = self.canvas.mpl_connect('key_press_event',self.plot_spectrum)

    def _disconnect(self):
        """ Disconnect click, click up (release click), and key press from events """
        self.canvas.mpl_disconnect(self.clickid)
        self.canvas.mpl_disconnect(self.clickupid)
        self.canvas.mpl_disconnect(self.keyid)

    def makeplane(self, estimator=np.mean):
        """

        *estimator* [ function | 'max' | 'int' | FITS filename | integer ]
            A non-pythonic, non-duck-typed variable.  If it's a function, apply that function
            along the cube's spectral axis to obtain an estimate (e.g., mean, min, max, etc.).
            'max' will do the same thing as passing np.max
            'int' will attempt to integrate the image (which is why I didn't duck-type)
            a .fits filename will be read using pyfits (so you can make your own cover figure)
            an integer will get the n'th slice in the parcube if it exists

        """
        # THIS IS A HACK!!!  isinstance(a function, function) must be a thing...
        FUNCTION = type(np.max)

        # estimator is NOT duck-typed
        if type(estimator) is FUNCTION:
            self.plane = estimator(self.Cube.cube,axis=0)
        elif type(estimator) is str:
            if estimator == 'max':
                self.plane = self.Cube.cube.max(axis=0)
            elif estimator == 'int':
                dx = np.abs(self.Cube.xarr[1:] - self.Cube.xarr[:-1])
                dx = np.concatenate([dx,[dx[-1]]])
                self.plane = (self.Cube.cube * dx[:,np.newaxis,np.newaxis]).sum(axis=0)
            elif estimator[-5:] == ".fits":
                self.plane = pyfits.getdata(estimator)
        elif type(estimator) is int:
            if hasattr(self.Cube,'parcube'):
                self.plane = self.Cube.parcube[estimator,:,:]

        if self.plane is None:
            raise ValueError("Invalid estimator %s" % (str(estimator)))

    def click(self,event):
        """
        Record location of downclick
        """
        if event.inaxes:
            self._clickX = event.xdata
            self._clickY = event.ydata

    def plot_spectrum(self, event, plot_fit=True):
        """
        Connects map cube to Spectrum...
        """
        self.event = event
        if event.inaxes:
            clickX = event.xdata
            clickY = event.ydata
        
            # grab toolbar info so that we don't do anything if a tool is selected
            tb = self.canvas.toolbar
            if tb.mode != '':
                return
            elif event.key is not None:
                if event.key == 'c':
                    self._center = (clickX,clickY)
                    self._remove_circle()
                    self._add_click_mark(clickX,clickY,clear=True)
                elif event.key == 'r':
                    x,y = self._center
                    self._add_circle(x,y,clickX,clickY)
                    self.circle(x,y,clickX,clickY)
            elif (hasattr(event,'button') and event.button in (1,2) 
                    and not (self._clickX == clickX and self._clickY == clickY)):
                if event.button == 1:
                    self._remove_circle()
                    clear=True
                    color = 'k'
                    linestyle = 'steps-mid'
                else:
                    color = self.overplot.colorcycle.next()
                    linestyle = self.overplot_linestyle
                    clear=False
                rad = ( (self._clickX-clickX)**2 + (self._clickY-clickY)**2 )**0.5
                print "Plotting circle from point %i,%i to %i,%i (r=%f)" % (self._clickX,self._clickY,clickX,clickY,rad)
                self._add_circle(self._clickX,self._clickY,clickX,clickY)
                self.circle(self._clickX,self._clickY,clickX,clickY,clear=clear,linestyle=linestyle,color=color)
            elif hasattr(event,'button') and event.button is not None:
                if event.button==1:
                    print "Plotting spectrum from point %i,%i" % (clickX,clickY)
                    self._remove_circle()
                    self._add_click_mark(clickX,clickY,clear=True)
                    self.Cube.plot_spectrum(clickX,clickY,clear=True)
                    if plot_fit: self.Cube.plot_fit(clickX, clickY, silent=True)
                elif event.button==2:
                    print "OverPlotting spectrum from point %i,%i" % (clickX,clickY)
                    color=self.overplot_colorcycle.next()
                    self._add_click_mark(clickX,clickY,clear=False, color=color)
                    self.Cube.plot_spectrum(clickX,clickY,clear=False, color=color, linestyle=self.overplot_linestyle)
                elif event.button==3:
                    print "Disconnecting GAIA-like tool"
                    self._disconnect()
            else:
                print "Call failed for some reason: "
                print "event: ",event
        else:
            warn("Click outside of axes")

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
            self._circles.append( matplotlib.patches.Circle([x,y],radius=r,**kwargs) )
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
                self.axis.patches.remove(circle)
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
