"""
=======
Plotter
=======

.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>
"""
from __future__ import print_function
import matplotlib
import matplotlib.figure
import numpy as np
import astropy.units as u
import copy
import inspect
from astropy import log

# this mess is to handle a nested hell of different versions of matplotlib
# (>=1.3 has BoundMethodProxy somewhere, >=3 gets rid of it) and python
# (python >=3.4 has WeakMethod, earlier versions don't)
try:
    from matplotlib.cbook import BoundMethodProxy
except ImportError:
    try:
        from matplotlib.cbook import _BoundMethodProxy as BoundMethodProxy
    except ImportError:
        try:
            from matplotlib.cbook import WeakMethod
        except ImportError:
            try:
                from weakref import WeakMethod
            except ImportError:
                try:
                    from weakrefmethod import WeakMethod
                except ImportError:
                    raise ImportError("Could not import WeakMethod from "
                                      "anywhere.  Try installing the "
                                      "weakrefmethod package or use a more "
                                      "recent version of python or matplotlib")

        class BoundMethodProxy(WeakMethod):
            @property
            def func(self):
                return self()

from . import widgets
from ..specwarnings import warn

interactive_help_message = """
Interactive key commands for plotter.  An additional help message may appear if
you have initiated the fitter.
'?' - bring up this message
'f' - initiate the /f/itter
'b' - initiate the /b/aseliner
'B' - initiate the /b/aseliner (reset the selection too)
'r' - re-attach matplotlib keys
'R' - redraw the plot cleanly
'i' : individual components / show each fitted component
"""

xlabel_table = {'speed': 'Velocity'}

class Plotter(object):
    """
    Class to plot a spectrum
    """


    def __init__(self, Spectrum, autorefresh=True, title="", xlabel=None,
                 silent=True, plotscale=1.0, **kwargs):

        import matplotlib.pyplot
        self._pyplot = matplotlib.pyplot

        self.figure = None
        self.axis = None
        self.Spectrum = Spectrum
        # plot parameters
        self.offset = 0.0 # vertical offset
        self.autorefresh = autorefresh
        self.xlabel = xlabel
        self.title = title
        self.errorplot = None
        self.plotkwargs = kwargs
        self._xlim = [None,None]
        self._ylim = [None,None]
        self.debug = False

        self.keyclick = None
        self.silent = silent
        self.plotscale = plotscale

        self._xclick1 = None
        self._xclick2 = None

        self.automake_fitter_tool = False

        self._active_gui = None

    @property
    def _xunit(self):
        return self.Spectrum.xarr.unit

    def _get_prop(xy, minmax):
        def getprop(self):
            if xy == 'x':
                if minmax == 'min':
                    if self._xlim[0] is not None and self._xunit:
                        try:
                            self._xlim[0]._unit = self._xunit
                        except AttributeError:
                            self._xlim[0] = u.Quantity(self._xlim[0], self._xunit)
                    return self._xlim[0]
                elif minmax == 'max':
                    if self._xlim[1] is not None and self._xunit:
                        try:
                            self._xlim[1]._unit = self._xunit
                        except AttributeError:
                            self._xlim[1] = u.Quantity(self._xlim[1], self._xunit)
                    return self._xlim[1]
            elif xy == 'y':
                if minmax == 'min':
                    return self._ylim[0]
                elif minmax == 'max':
                    return self._ylim[1]
        return getprop

    def _set_prop(xy, minmax):
        def setprop(self, value):
            if self.debug:
                frm = inspect.stack()
                print(frm[1],"Setting %s%s to %s" % (xy,minmax,value))
            if xy == 'x':
                if minmax == 'min':
                    self._xlim[0] = value
                elif minmax == 'max':
                    self._xlim[1] = value
            elif xy == 'y':
                if minmax == 'min':
                    self._ylim[0] = value
                elif minmax == 'max':
                    self._ylim[1] = value
        return setprop

    xmin = property(fget=_get_prop('x','min'),fset=_set_prop('x','min'))
    xmax = property(fget=_get_prop('x','max'),fset=_set_prop('x','max'))
    ymin = property(fget=_get_prop('y','min'),fset=_set_prop('y','min'))
    ymax = property(fget=_get_prop('y','max'),fset=_set_prop('y','max'))


    def _disconnect_matplotlib_keys(self):
        """
        Disconnected the matplotlib key-press callbacks
        """
        if self.figure is not None:
            cbs = self.figure.canvas.callbacks.callbacks
            # this may cause problems since the dict of key press events is a
            # dict, i.e. not ordered, and we want to pop the first one...
            mpl_keypress_handler = self.figure.canvas.manager.key_press_handler_id
            try:
                self._mpl_key_callbacks = {mpl_keypress_handler:
                                           cbs['key_press_event'].pop(mpl_keypress_handler)}
            except KeyError:
                bmp = BoundMethodProxy(self.figure.canvas.manager.key_press)
                self._mpl_key_callbacks = {mpl_keypress_handler:
                                           bmp}

    def _reconnect_matplotlib_keys(self):
        """
        Reconnect the previously disconnected matplotlib keys
        """
        if self.figure is not None and hasattr(self,'_mpl_key_callbacks'):
            self.figure.canvas.callbacks.callbacks['key_press_event'].update(self._mpl_key_callbacks)
        elif self.figure is not None:
            mpl_keypress_handler = self.figure.canvas.manager.key_press_handler_id
            try:
                bmp = BoundMethodProxy(self.figure.canvas.manager.key_press)
                self.figure.canvas.callbacks.callbacks['key_press_event'].update({mpl_keypress_handler:
                                                                                  bmp})
            except AttributeError as ex:
                print(f"Error {ex} was raised when trying to connect the key_press handler.  "
                      "Please file an issue on github.  You may try a different matplotlib backend "
                      "as a temporary workaround")

    def __call__(self, figure=None, axis=None, clear=True, autorefresh=None,
                 plotscale=1.0, override_plotkwargs=False, **kwargs):
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
            self.figure = self._pyplot.figure(figure)
            self.axis = self.figure.gca()
        elif self.figure is None:
            if isinstance(axis,matplotlib.axes.Axes):
                self.axis = axis
                self.figure = axis.figure
            else:
                self.figure = self._pyplot.figure()

        if hasattr(self.figure, 'number') and not self._pyplot.fignum_exists(self.figure.number):
            self.figure = self._pyplot.figure(self.figure.number)

        # always re-connect the interactive keys to avoid frustration...
        self._mpl_reconnect()

        if axis is not None:
            #self._mpl_disconnect()
            self.axis = axis
            self.figure = axis.figure
            #self._mpl_connect()
        elif len(self.figure.axes) > 0 and self.axis is None:
            self.axis = self.figure.axes[0] # default to first axis
        elif self.axis is None:
            self.axis = self.figure.gca()

        # A check to deal with issue #117: if you close the figure, the axis
        # still exists, but it cannot be reattached to a figure
        if (hasattr(self.axis.get_figure(), 'number') and
            not (self.axis.get_figure() is self._pyplot.figure(self.axis.get_figure().number))):
            self.axis = self.figure.gca()

        if self.axis is not None and self.axis not in self.figure.axes:
            # if you've cleared the axis, but the figure is still open, you
            # need a new axis
            self.figure.add_axes(self.axis)


        if clear and self.axis is not None:
            self.axis.clear()
            # Need to empty the stored model plots
            if hasattr(self.Spectrum, 'fitter'):
                self.Spectrum.fitter.clear()

        if autorefresh is not None:
            self.autorefresh = autorefresh

        self.plotscale = plotscale

        if self.plotkwargs and not override_plotkwargs:
            self.plotkwargs.update(kwargs)
        else:
            self.plotkwargs = kwargs

        self.plot(**kwargs)

    def _mpl_connect(self):
        if self.keyclick is None:
            self.keyclick = self.figure.canvas.mpl_connect('key_press_event',self.parse_keys)

    def _mpl_disconnect(self):
        self.figure.canvas.mpl_disconnect(self.keyclick)
        self.keyclick = None

    def disconnect(self):
        """
        Disconnect the matplotlib interactivity of this pyspeckit plotter.
        """
        self._mpl_disconnect()

    def connect(self):
        """
        Connect to the matplotlib key-parsing interactivity
        """
        self._mpl_connect()

    def _mpl_reconnect(self):
        self._mpl_disconnect()
        self._mpl_connect()
        # disable fullscreen & grid
        self._pyplot.rcParams['keymap.fullscreen'] = 'ctrl+f'
        self._pyplot.rcParams['keymap.grid'] = 'ctrl+g'

    def plot(self, offset=0.0, xoffset=0.0, color='k', drawstyle='steps-mid',
             linewidth=0.5, errstyle=None, erralpha=0.2, errcolor=None,
             silent=None, reset=True, refresh=True, use_window_limits=None,
             useOffset=False, **kwargs):
        """
        Plot the spectrum!

        Tries to automatically find a reasonable plotting range if one is not
        set.

        Parameters
        ----------
        offset : float
            vertical offset to add to the spectrum before plotting.  Useful if
            you want to overlay multiple spectra on a single plot
        xoffset: float
            An x-axis shift.  I don't know why you'd want this...
        color : str
            default to plotting spectrum in black
        drawstyle : 'steps-mid' or str
            'steps-mid' for histogram-style plotting.  See matplotlib's plot
            for more information
        linewidth : float
            Line width in pixels.  Narrow lines are helpful when histo-plotting
        errstyle : 'fill', 'bars', or None
            can be "fill", which draws partially transparent boxes around the
            data to show the error region, or "bars" which draws standard
            errorbars.  ``None`` will display no errorbars
        useOffset : bool
            Use offset-style X/Y coordinates (e.g., 1 + 1.483e10)?  Defaults to
            False because these are usually quite annoying.
        xmin/xmax/ymin/ymax : float
            override defaults for plot range.  Once set, these parameters are
            sticky (i.e., replotting will use the same ranges).  Passed to
            `reset_limits`
        reset_[xy]limits : bool
            Reset the limits to "sensible defaults".  Passed to `reset_limits`
        ypeakscale : float
            Scale up the Y maximum value.  Useful to keep the annotations away
            from the data.  Passed to `reset_limits`
        reset : bool
            Reset the x/y axis limits? If set, `reset_limits` will be called.
        """

        if self.axis is None:
            raise Exception("You must call the Plotter class to initiate the canvas before plotting.")

        self.offset = offset

        # there is a bug where this only seems to update the second time it is called
        self.label(**kwargs)
        self.label(**kwargs)
        for arg in ['title','xlabel','ylabel']:
            if arg in kwargs:
                kwargs.pop(arg)

        reset_kwargs = {}
        for arg in ['xmin', 'xmax', 'ymin', 'ymax', 'reset_xlimits',
                    'reset_ylimits', 'ypeakscale']:
            if arg in kwargs:
                reset_kwargs[arg] = kwargs.pop(arg)

        if (use_window_limits is None and any(k in reset_kwargs for k in
                                              ('xmin','xmax','reset_xlimits'))):
            use_window_limits = False

        if use_window_limits:
            self._stash_window_limits()

        # for filled errorbars, order matters.
        inds = np.argsort(self.Spectrum.xarr)

        if errstyle is not None:
            if errcolor is None:
                errcolor = color
            if errstyle == 'fill':
                self.errorplot = [self.axis.fill_between(steppify(self.Spectrum.xarr.value[inds]+xoffset, isX=True),
                                  steppify((self.Spectrum.data*self.plotscale+self.offset-self.Spectrum.error*self.plotscale)[inds]),
                                  steppify((self.Spectrum.data*self.plotscale+self.offset+self.Spectrum.error*self.plotscale)[inds]),
                                  facecolor=errcolor, edgecolor=errcolor, alpha=erralpha, **kwargs)]
            elif errstyle == 'bars':
                self.errorplot = self.axis.errorbar(self.Spectrum.xarr[inds].value+xoffset,
                                                    self.Spectrum.data[inds]*self.plotscale+self.offset,
                                                    yerr=self.Spectrum.error[inds]*self.plotscale,
                                                    ecolor=errcolor, fmt='none',
                                                    **kwargs)

        self._spectrumplot = self.axis.plot(self.Spectrum.xarr.value[inds]+xoffset,
                                            self.Spectrum.data[inds]*self.plotscale+self.offset,
                                            color=color,
                                            drawstyle=drawstyle,
                                            linewidth=linewidth, **kwargs)

        self.axis.ticklabel_format(useOffset=useOffset)

        if use_window_limits:
            self._reset_to_stashed_limits()

        if silent is not None:
            self.silent = silent

        if reset:
            self.reset_limits(use_window_limits=use_window_limits, **reset_kwargs)

        if self.autorefresh and refresh:
            self.refresh()

        # Maybe it's OK to call 'plot' when there is an active gui tool
        # (e.g., baseline or specfit)?
        #if self._active_gui:
        #    self._active_gui = None
        #    warn("An active GUI was found while initializing the "
        #         "plot.  This is somewhat dangerous and may result "
        #         "in broken interactivity.")


    def _stash_window_limits(self):
        self._window_limits = self.axis.get_xlim(),self.axis.get_ylim()
        if self.debug:
            print("Stashed window limits: ",self._window_limits)

    def _reset_to_stashed_limits(self):
        self.axis.set_xlim(*self._window_limits[0])
        self.axis.set_ylim(*self._window_limits[1])
        self.xmin,self.xmax = self._window_limits[0]
        self.ymin,self.ymax = self._window_limits[1]
        if self.debug:
            print("Recovered window limits: ",self._window_limits)

    def reset_limits(self, xmin=None, xmax=None, ymin=None, ymax=None,
                     reset_xlimits=True, reset_ylimits=True, ypeakscale=1.2,
                     silent=None, use_window_limits=False, **kwargs):
        """
        Automatically or manually reset the plot limits
        """
        # if not use_window_limits: use_window_limits = False
        if self.debug:
            frame = inspect.currentframe()
            args, _, _, values = inspect.getargvalues(frame)
            print(zip(args,values))

        if use_window_limits:
            # this means DO NOT reset!
            # it simply sets self.[xy][min/max] = current value
            self.set_limits_from_visible_window()
        else:
            if silent is not None:
                self.silent = silent

            # if self.xmin and self.xmax:
            if (reset_xlimits or self.Spectrum.xarr.min().value < self.xmin or self.Spectrum.xarr.max().value > self.xmax):
                if not self.silent:
                    warn("Resetting X-axis min/max because the plot is out of bounds.")
                self.xmin = None
                self.xmax = None
            if xmin is not None:
                self.xmin = u.Quantity(xmin, self._xunit)
            elif self.xmin is None:
                self.xmin = u.Quantity(self.Spectrum.xarr.min().value, self._xunit)
            if xmax is not None:
                self.xmax = u.Quantity(xmax, self._xunit)
            elif self.xmax is None:
                self.xmax = u.Quantity(self.Spectrum.xarr.max().value, self._xunit)

            xpixmin = np.argmin(np.abs(self.Spectrum.xarr.value-self.xmin.value))
            xpixmax = np.argmin(np.abs(self.Spectrum.xarr.value-self.xmax.value))
            if xpixmin>xpixmax:
                xpixmin,xpixmax = xpixmax,xpixmin
            elif xpixmin == xpixmax:
                if reset_xlimits:
                    raise Exception("Infinite recursion error.  Maybe there are no valid data?")
                if not self.silent:
                    warn("ERROR: the X axis limits specified were invalid.  Resetting.")
                self.reset_limits(reset_xlimits=True, ymin=ymin, ymax=ymax,
                                  reset_ylimits=reset_ylimits,
                                  ypeakscale=ypeakscale, **kwargs)
                return

            if self.ymin is not None and self.ymax is not None:
                # this is utter nonsense....
                if (np.nanmax(self.Spectrum.data) < self.ymin or np.nanmin(self.Spectrum.data) > self.ymax
                        or reset_ylimits):
                    if not self.silent and not reset_ylimits:
                        warn("Resetting Y-axis min/max because the plot is out of bounds.")
                    self.ymin = None
                    self.ymax = None

            if ymin is not None:
                self.ymin = ymin
            elif self.ymin is None:
                yminval = np.nanmin(self.Spectrum.data[xpixmin:xpixmax])
                # Increase the range fractionally.  This means dividing a positive #, multiplying a negative #
                if yminval < 0:
                    self.ymin = float(yminval)*float(ypeakscale)
                else:
                    self.ymin = float(yminval)/float(ypeakscale)

            if ymax is not None:
                self.ymax = ymax
            elif self.ymax is None:
                ymaxval = (np.nanmax(self.Spectrum.data[xpixmin:xpixmax])-self.ymin)
                if ymaxval > 0:
                    self.ymax = float(ymaxval) * float(ypeakscale) + self.ymin
                else:
                    self.ymax = float(ymaxval) / float(ypeakscale) + self.ymin

            self.ymin += self.offset
            self.ymax += self.offset

        self.axis.set_xlim(self.xmin.value if hasattr(self.xmin, 'value') else self.xmin,
                           self.xmax.value if hasattr(self.xmax, 'value') else self.xmax)
        self.axis.set_ylim(self.ymin, self.ymax)


    def label(self, title=None, xlabel=None, ylabel=None, verbose_label=False,
              **kwargs):
        """
        Label the plot, with an attempt to parse standard units into nice latex labels

        Parameters
        ----------
        title : str
        xlabel : str
        ylabel : str
        verbose_label: bool
        """

        if title is not None:
            self.title = title
        elif hasattr(self.Spectrum,'specname'):
            self.title = self.Spectrum.specname
        if self.title != "":
            self.axis.set_title(self.title)

        if xlabel is not None:
            log.debug("setting xlabel={0}".format(xlabel))
            self.xlabel = xlabel
        elif self._xunit:
            try:
                self.xlabel = xlabel_table[str(self._xunit.physical_type).lower()]
            except KeyError:
                self.xlabel = str(self._xunit.physical_type)
            # WAS: self.xlabel += " ("+u.Unit(self._xunit).to_string()+")"
            self.xlabel += " ({0})".format(self._xunit.to_string())
            log.debug("xunit is {1}. set xlabel={0}".format(self.xlabel,
                                                            self._xunit))

            if verbose_label:
                self.xlabel = "%s %s" % (str(self.Spectrum.xarr.velocity_convention),
                                         self.xlabel)
        else:
            log.warn("Plotter: xlabel was not set")

        if self.xlabel is not None:
            self.axis.set_xlabel(self.xlabel)

        if ylabel is not None:
            self.axis.set_ylabel(ylabel)
        elif self.Spectrum.unit in ['Ta*','Tastar']:
            self.axis.set_ylabel("$T_A^*$ (K)")
        elif self.Spectrum.unit in ['K']:
            self.axis.set_ylabel("Brightness Temperature $T$ (K)")
        elif self.Spectrum.unit == 'mJy':
            self.axis.set_ylabel("$S_\\nu$ (mJy)")
        elif self.Spectrum.unit == 'Jy':
            self.axis.set_ylabel("$S_\\nu$ (Jy)")
        else:
            if isinstance(self.Spectrum.unit, str) and "$" in self.Spectrum.unit:
                # assume LaTeX already
                self.axis.set_ylabel(self.Spectrum.unit)
            elif isinstance(self.Spectrum.unit, str):
                self.axis.set_ylabel(self.Spectrum.unit)
            else:
                label_units = self.Spectrum.unit.to_string(format='latex')
                if 'mathring{A}' in label_units:
                    label_units = label_units.replace('\\mathring{A}', 'A')
                if '\\overset' in label_units:
                    label_units = label_units.replace('\\overset', '^')
                self.axis.set_ylabel(label_units)

    @property
    def ylabel(self):
        return self.axis.get_ylabel()

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
                print(interactive_help_message)
            elif event.key == 'f':
                print("\n\nFitter initiated from the interactive plotter.")
                # extra optional text:
                #  Matplotlib shortcut keys ('g','l','p',etc.) are disabled.  Re-enable with 'r'"
                if self._active_gui == self.Spectrum.specfit and self._active_gui._check_connections(verbose=False):
                    print("Fitter is already active.  Use 'q' to quit the fitter.")
                elif self._active_gui == self.Spectrum.specfit and not self._active_gui._check_connections(verbose=False):
                    # forcibly clear connections
                    self._active_gui.clear_all_connections()
                    # the 'clear_all_connections' code *explicitly* makes the
                    # following line correct, except in the case that there is
                    # no canvas...
                    assert self._active_gui is None
                    self.activate_interactive_fitter()
                else:
                    self.activate_interactive_fitter()

                assert self._active_gui == self.Spectrum.specfit
                assert self._active_gui._check_connections(verbose=False)

                if not hasattr(self,'FitterTool') and self.automake_fitter_tool:
                    self.FitterTool = widgets.FitterTools(self.Spectrum.specfit, self.figure)
                elif hasattr(self,'FitterTool') and self.FitterTool.toolfig.number not in self._pyplot.get_fignums():
                    self.FitterTool = widgets.FitterTools(self.Spectrum.specfit, self.figure)
            elif event.key is not None and event.key.lower() == 'b':
                if event.key == 'b':
                    print("\n\nBaseline initiated from the interactive plotter")
                elif event.key == 'B':
                    print("\n\nBaseline initiated from the interactive plotter (with reset)")
                print("Matplotlib shortcut keys ('g','l','p',etc.) are disabled.  Re-enable with 'r'")
                self.activate_interactive_baseline_fitter(reset_selection=(event.key=='B'))

                if not hasattr(self,'FitterTool') and self.automake_fitter_tool:
                    self.FitterTool = widgets.FitterTools(self.Spectrum.specfit, self.figure)
                elif hasattr(self,'FitterTool') and self.FitterTool.toolfig.number not in self._pyplot.get_fignums():
                    self.FitterTool = widgets.FitterTools(self.Spectrum.specfit, self.figure)
            elif event.key == 'r':
                # print("\n\nReconnected matplotlib shortcut keys.")
                self._reconnect_matplotlib_keys()
            elif event.key == 'R':
                self()
            elif event.key == 'i':
                self.Spectrum.specfit.plot_fit(show_components=True)

    def get_two_clicks(self,event):

        if self._xclick1 is None:
            self._xclick1 = event.xdata
        elif self._xclick2 is None:
            self._xclick2 = event.xdata

    def set_limits_from_visible_window(self, debug=False):
        """ Hopefully self-descriptive: set the x and y limits from the
        currently visible window (use this if you use the pan/zoom tools or
        manually change the limits) """
        if debug:
            print("Changing x limits from {},{} to {},{}".format(self.xmin,self.xmax,self.axis.get_xlim()[0],self.axis.get_xlim()[1]))
            print("Changing y limits from {},{} to {},{}".format(self.ymin,self.ymax,self.axis.get_ylim()[0],self.axis.get_ylim()[1]))
        self.xmin, self.xmax = self.axis.get_xlim()
        self.ymin, self.ymax = self.axis.get_ylim()
        if debug:
            print("New x limits {},{} == {},{}".format(self.xmin,self.xmax,self.axis.get_xlim()[0],self.axis.get_xlim()[1]))
            print("New y limits {},{} == {},{}".format(self.ymin,self.ymax,self.axis.get_ylim()[0],self.axis.get_ylim()[1]))

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

    def line_ids(self, line_names, line_xvals, xval_units=None, auto_yloc=True,
                 velocity_offset=None, velocity_convention='radio',
                 auto_yloc_fraction=0.9,  **kwargs):
        """
        Add line ID labels to a plot using lineid_plot
        http://oneau.wordpress.com/2011/10/01/line-id-plot/
        https://github.com/phn/lineid_plot
        http://packages.python.org/lineid_plot/

        Parameters
        ----------
        line_names : list
            A list of strings to label the specified x-axis values
        line_xvals : list
            List of x-axis values (e.g., wavelengths) at which to label the lines.
            Can be a list of quantities.
        xval_units : string
            The unit of the line_xvals if they are not given as quantities
        velocity_offset : quantity
            A velocity offset to apply to the inputs if they are in frequency
            or wavelength units
        velocity_convention : 'radio' or 'optical' or 'doppler'
            Used if the velocity offset is given
        auto_yloc : bool
            If set, overrides box_loc and arrow_tip (the vertical position of
            the lineid labels) in kwargs to be `auto_yloc_fraction` of the plot
            range
        auto_yloc_fraction: float in range [0,1]
            The fraction of the plot (vertically) at which to place labels

        Examples
        --------
        >>> import numpy as np
        >>> import pyspeckit
        >>> sp = pyspeckit.Spectrum(
                xarr=pyspeckit.units.SpectroscopicAxis(np.linspace(-50,50,101),
                    unit='km/s', refX=6562.8, refX_unit='angstrom'),
                data=np.random.randn(101), error=np.ones(101))
        >>> sp.plotter()
        >>> sp.plotter.line_ids(['H$\\alpha$'],[6562.8],xval_units='angstrom')
        """
        import lineid_plot

        if velocity_offset is not None:
            assert velocity_offset.unit.is_equivalent(u.km/u.s)

        doppler = getattr(u, 'doppler_{0}'.format(velocity_convention))
        if self.Spectrum.xarr.refX is not None:
            equivalency = doppler(self.Spectrum.xarr.refX)
        else:
            equivalency = doppler(self.Spectrum.xarr.as_unit(u.GHz)[0])

        xvals = []
        linenames_toplot = []
        for xv,ln in zip(line_xvals, line_names):
            if hasattr(xv, 'unit'):
                pass
            else:
                xv = u.Quantity(xv, xval_units)

            xv = xv.to(u.km/u.s,
                       equivalencies=equivalency)
            if velocity_offset is not None:
                xv = xv + velocity_offset
            xv = xv.to(self.Spectrum.xarr.unit, equivalencies=equivalency)

            if self.Spectrum.xarr.in_range(xv):
                xvals.append(xv.value)
                linenames_toplot.append(ln)

        if len(xvals) != len(line_xvals):
            log.warn("Skipped {0} out-of-bounds lines when plotting line IDs."
                     .format(len(line_xvals)-len(xvals)))

        if auto_yloc:
            yr = self.axis.get_ylim()
            kwargs['box_loc'] = (yr[1]-yr[0])*auto_yloc_fraction + yr[0]
            kwargs['arrow_tip'] = (yr[1]-yr[0])*(auto_yloc_fraction*0.9) + yr[0]

        lineid_plot.plot_line_ids(self.Spectrum.xarr,
                                  self.Spectrum.data,
                                  xvals,
                                  linenames_toplot,
                                  ax=self.axis,
                                  **kwargs)

    def line_ids_from_measurements(self, auto_yloc=True,
                                   auto_yloc_fraction=0.9, **kwargs):
        """
        Add line ID labels to a plot using lineid_plot
        http://oneau.wordpress.com/2011/10/01/line-id-plot/
        https://github.com/phn/lineid_plot
        http://packages.python.org/lineid_plot/

        Parameters
        ----------
        auto_yloc : bool
            If set, overrides box_loc and arrow_tip (the vertical position of
            the lineid labels) in kwargs to be `auto_yloc_fraction` of the plot
            range
        auto_yloc_fraction: float in range [0,1]
            The fraction of the plot (vertically) at which to place labels

        Examples
        --------
        >>> import numpy as np
        >>> import pyspeckit
        >>> sp = pyspeckit.Spectrum(
                xarr=pyspeckit.units.SpectroscopicAxis(np.linspace(-50,50,101),
                    units='km/s', refX=6562.8, refX_unit='angstroms'),
                data=np.random.randn(101), error=np.ones(101))
        >>> sp.plotter()
        >>> sp.specfit(multifit=None, fittype='gaussian', guesses=[1,0,1]) # fitting noise....
        >>> sp.measure()
        >>> sp.plotter.line_ids_from_measurements()
        """
        import lineid_plot

        if hasattr(self.Spectrum,'measurements'):
            measurements = self.Spectrum.measurements

            if auto_yloc:
                yr = self.axis.get_ylim()
                kwargs['box_loc'] = (yr[1]-yr[0])*auto_yloc_fraction + yr[0]
                kwargs['arrow_tip'] = (yr[1]-yr[0])*(auto_yloc_fraction*0.9) + yr[0]

            lineid_plot.plot_line_ids(self.Spectrum.xarr, self.Spectrum.data,
                                      [v['pos'] for v in
                                       measurements.lines.values()],
                                      measurements.lines.keys(), ax=self.axis,
                                      **kwargs)
        else:
            warn("Cannot add line IDs from measurements unless measurements have been made!")

    def activate_interactive_fitter(self):
        """
        Attempt to activate the interactive fitter
        """
        if self._active_gui is not None:
            # This should not be reachable.  Clearing connections is the
            # "right" behavior if this becomes reachable, but I'd rather raise
            # an exception because I don't want to get here ever
            self._active_gui.clear_all_connections()
            raise ValueError("GUI was active when 'f' key pressed")

        self._activate_interactive(self.Spectrum.specfit, interactive=True)

    def activate_interactive_baseline_fitter(self, **kwargs):
        """
        Attempt to activate the interactive baseline fitter
        """
        if self._active_gui is not None:
            # This should not be reachable.  Clearing connections is the
            # "right" behavior if this becomes reachable, but I'd rather raise
            # an exception because I don't want to get here ever
            gui_was = self._active_gui
            self._active_gui.clear_all_connections()
            raise ValueError("GUI {0} was active when 'b' key pressed"
                             .format(gui_was))

        self._activate_interactive(self.Spectrum.baseline, interactive=True,
                                   **kwargs)

    def _activate_interactive(self, object_to_activate, **kwargs):
        self._disconnect_matplotlib_keys()

        self._active_gui = object_to_activate

        # activating the gui calls clear_all_connections, which disconnects the
        # gui
        try:
            self._active_gui(**kwargs)
            self._active_gui = object_to_activate
            assert self._active_gui is not None
        except Exception as ex:
            self._active_gui = None
            raise ex

def parse_units(labelstring):
    import re
    labelstring = re.sub("um","$\\mu$m",labelstring)
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
        newarr = np.array(list(zip(arr[:-1]-interval,arr[:-1]+interval))).ravel()
        newarr = np.concatenate([newarr,2*[newarr[-1]+interval[-1]]])
    else:
        newarr = np.array(list(zip(arr,arr))).ravel()
    return newarr
