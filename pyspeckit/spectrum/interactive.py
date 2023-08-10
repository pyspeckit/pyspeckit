"""
==================
Window Interaction
==================

A general module for selecting regions and inputting guesses via the
interactive window.
"""
from __future__ import print_function
import numpy
from . import units
from astropy import log
from six.moves import xrange
from six import iteritems

class Interactive(object):

    def __init__(self, Spectrum, guesses=None,
                 interactive_help_message="Replace this message"):
        """
        Declare interactive variables.

        Must have a parent Spectrum class

        **Must declare button2action and button3action**
        """
        self.Spectrum = Spectrum

        self.interactive_help_message = interactive_help_message
        # includemask should not be a masked array even if data is
        # masked arrays are apparently full of bugs...
        self.includemask = numpy.ones(self.Spectrum.data.size, dtype='bool')
        self.xclicks = []
        self.yclicks = []
        self.event_history = []
        self.guesses = guesses

        # Click counters
        self.nclicks_b1 = 0 # button 1
        self.nclicks_b2 = 0 # button 2

        # Temporary storage (for left, right clicking)
        self._xclick1 = None
        self._xclick2 = None

        # Set min/max of range
        self.xmin = 0
        self.xmax = self.Spectrum.xarr.shape[0]

        # Init button 1/2 plots
        self.button1plot = []
        self.button2plot = []

        self.use_window_limits = None

        # initialization: Glue can activate fitters without start_interactive,
        # so these need to be declared
        self.click = None
        self.keyclick = None

        self._debug = False

    @property
    def xmin(self):
        return self._xmin

    @xmin.setter
    def xmin(self, value):
        self._xmin = int(value)

    @property
    def xmax(self):
        return self._xmax

    @xmax.setter
    def xmax(self, value):
        self._xmax = int(value)

    def event_manager(self, event, force_over_toolbar=False, debug=False):
        """
        Decide what to do given input (click, keypress, etc.)
        """
        if hasattr(self.Spectrum.plotter.figure.canvas.manager, 'toolbar'):
            toolbar = self.Spectrum.plotter.figure.canvas.manager.toolbar
            toolmode = toolbar.mode
        else:
            # If interactivity isn't possible, we don't really care what tool is 'active'
            toolmode = ''
        self.event_history.append(event)

        #DEBUG print("toolmode = {0} force_over_toolbar={1}".format(toolmode, force_over_toolbar))
        if (toolmode == '' or force_over_toolbar) and self.Spectrum.plotter.axis in event.canvas.figure.axes:
            if hasattr(event,'button'):
                button = event.button
            elif hasattr(event,'key'):
                button = event.key

            #DEBUG print("Event: {0}".format(event))
            if event.xdata is None or event.ydata is None:
                return

            if debug or self._debug:
                log.debug("button: {0} x,y: {1},{2} "
                          " nclicks 1: {3:f}  2: {4:f}".format
                          (button, event.xdata, event.ydata, self.nclicks_b1,
                           self.nclicks_b2,))

            if button in ('p','P','1',1,'i','a'): # p for... parea?  a for area.  i for include
                # button one is always region selection
                #DEBUG print("Button is {0}".format(button))
                self._selectregion_interactive(event,debug=debug)
            elif button in ('c','C'):
                self.clear_highlights()
                self.clear_all_connections()
                self.Spectrum.plotter()
            elif button in ('e','x','E','X'): # e for exclude, x for x-clude
                # exclude/delete/remove
                self._selectregion_interactive(event, mark_include=False, debug=debug)
            elif button in ('m','M','2',2): # m for mark
                if debug or self._debug:
                    log.debug("Button 2 action")
                self.button2action(event,debug=debug)
            elif button in ('d','D','3',3): # d for done
                if debug or self._debug:
                    log.debug("Button 3 action")
                self.button3action(event,debug=debug)
            elif button in ('?'):
                # print statement: we really want this to go to the terminal
                print(self.interactive_help_message)
            elif hasattr(self,'Registry') and button in self.Registry.fitkeys:
                fittername = self.Registry.fitkeys[button]
                if fittername in self.Registry.multifitters:
                    self.fitter = self.Registry.multifitters[fittername]
                    self.fittype = fittername
                    print("Selected multi-fitter %s" % fittername)
                else:
                    print("ERROR: Did not find fitter %s" % fittername)
            if self.Spectrum.plotter.autorefresh: self.Spectrum.plotter.refresh()
        elif debug or self._debug:
            print("Button press not acknowledged.  event={0}, toolmode={1}".format(event,
                                                                                   toolmode))


    def _selectregion_interactive(self, event, mark_include=True, debug=False, **kwargs):
        """
        select regions for baseline fitting
        """

        xpix = self.Spectrum.xarr.x_to_pix(event.xdata)
        if self.xclicks == []:
            self._firstclick_selection(not mark_include)

        if self.nclicks_b1 == 0:
            self.nclicks_b1 = 1
            self._xclick1 = xpix
            self.xclicks.append(xpix)
            if debug or self._debug:
                print("Click 1: clickx=%i xmin=%i, xmax=%i" % (xpix,self.xmin,self.xmax))
        elif self.nclicks_b1 == 1:
            self._xclick2 = xpix
            self.nclicks_b1 = 0
            self.xclicks.append(xpix)
            # force click1 to be left (swap)
            if self._xclick1 > self._xclick2:
                self._xclick1,self._xclick2 = self._xclick2,self._xclick1

            # ensure that the fit/plot range is at least as large as the click range
            if self.xmin > self._xclick1: self.xmin = self._xclick1
            if self.xmax < self._xclick2: self.xmax = self._xclick2

            # change the includemask
            self.includemask[self._xclick1:self._xclick2] = mark_include

            if mark_include:
                self.highlight_fitregion(**kwargs)
            else:  # mark include=False -> mark_exclude=True
                for highlight_line in self.button1plot:
                    hlx,hly = highlight_line.get_data()
                    hide = ((hlx > self.Spectrum.xarr[self._xclick1]) *
                            (hlx < self.Spectrum.xarr[self._xclick2]))
                    hly[hide] = numpy.nan
                    highlight_line.set_ydata(hly)
                self.Spectrum.plotter.refresh()

            if debug or self._debug:
                print("Click 2: clickx=%i xmin=%i, xmax=%i" % (xpix,self.xmin,self.xmax))

        self._update_xminmax()

    def highlight_fitregion(self,  drawstyle='steps-mid', color=(0,0.8,0,0.5),
                            linewidth=2, alpha=0.5, clear_highlights=True,
                            **kwargs):
        """
        Re-highlight the fitted region

        kwargs are passed to `matplotlib.plot`
        """

        if clear_highlights:
            self.clear_highlights()

        bad = self.Spectrum.data*0
        bad[~self.includemask] = numpy.nan

        self.button1plot += self.Spectrum.plotter.axis.plot(
                self.Spectrum.xarr,
                # +bad adds nans to points that are not to be included
                self.Spectrum.data+self.Spectrum.plotter.offset+bad,
                drawstyle=drawstyle, color=color,
                linewidth=linewidth,
                alpha=alpha,
                **kwargs)
        self.Spectrum.plotter.refresh()

    def _firstclick_selection(self, include_all=False):
        """
        Initialize the include/exclude mask
        """

        self.Spectrum.plotter.axis.set_autoscale_on(False)
        if include_all:
            # default to including everything
            self.includemask = numpy.array(self.Spectrum.data, dtype='bool') + True
        else:
            # default to including nothing
            self.includemask = numpy.array(self.Spectrum.data, dtype='bool') * False

    def guesspeakwidth(self,event,debug=False,nwidths=1,**kwargs):
        """
        Interactively guess the peak height and width from user input

        Width is assumed to be half-width-half-max
        """
        modnum = 1+nwidths
        if debug or self._debug: print("nclicks: %i nwidths: %i modnum: %i" % (self.nclicks_b2,nwidths,modnum))
        if self.nclicks_b2 == 0:
            self.firstclick_guess()
        if self.nclicks_b2 % modnum == 0:
            # even clicks are peaks
            if self.Spectrum.baseline.subtracted:
                peakguess = event.ydata
            else:
                peakguess = event.ydata - self.Spectrum.baseline.basespec[self.Spectrum.xarr.x_to_pix(event.xdata)]
            self.guesses += [peakguess,event.xdata] + [1]*nwidths
            self.npeaks += 1
            self.nclicks_b2 += 1
            if debug or self._debug:
                print("Peak %i click %i at x,y %g,%g" % (self.npeaks,self.nclicks_b2,event.xdata,event.ydata))
            self.button2plot += [self.Spectrum.plotter.axis.scatter(event.xdata,event.ydata,marker='x',c='r')]
            #self.Spectrum.plotter.refresh() #plot(**self.Spectrum.plotter.plotkwargs)
        elif self.nclicks_b2 % modnum >= 1:
            # odd clicks are widths
            whichwidth = self.nclicks_b2 % modnum
            widthguess = (abs(event.xdata-self.guesses[-1-nwidths]) /
                          numpy.sqrt(2*numpy.log(2)))
            if numpy.isnan(widthguess) or widthguess <= 0:
                newwidthguess = numpy.abs(self.Spectrum.xarr.diff()).min().value
                if newwidthguess <= 0:
                    raise ValueError("A width guess could not be determined.")
                log.exception("Error: width guess was {0}.  It is being forced to {1}."
                              .format(widthguess, newwidthguess))
                widthguess = newwidthguess
            self.guesses[-whichwidth] = widthguess
            if debug or self._debug:
                print("Width %i whichwidth %i click %i at x,y %g,%g width: %g" % (self.npeaks,whichwidth,self.nclicks_b2,event.xdata,event.ydata,self.guesses[-whichwidth]))
            self.button2plot += self.Spectrum.plotter.axis.plot([event.xdata,
                                                                 2*self.guesses[-1-nwidths]-event.xdata],
                                                                [event.ydata]*2,
                                                                color='r')
            #self.Spectrum.plotter.refresh() #plot(**self.Spectrum.plotter.plotkwargs)
            if self.nclicks_b2 / (1+nwidths) > self.npeaks:
                print("There have been %i middle-clicks but there are only %i features" % (self.nclicks_b2,self.npeaks))
                self.npeaks += 1
            self.nclicks_b2 += 1
        else:
            raise ValueError("Bug in guesspeakwidth: somehow, the number of clicks doesn't make sense.")
        if debug or self._debug:
            print("Guesses: ",self.guesses)

    def firstclick_guess(self):
        """
        Initialize self.guesses
        """
        self.Spectrum.plotter.axis.set_autoscale_on(False)
        if self.guesses is None:
            self.guesses = []
        elif len(self.guesses) > 0:
            for ii in xrange(len(self.guesses)):
                self.guesses.pop()

    def clear_all_connections(self, debug=False):
        """
        Prevent overlapping interactive sessions
        """
        # this is really ugly, but needs to be done in order to prevent multiple overlapping calls...
        cids_to_remove = []
        if not hasattr(self.Spectrum.plotter.figure,'canvas'):
            # just quit out; saves a tab...
            if debug or self._debug:
                print("Didn't find a canvas, quitting.")
            # just in case?  This should be *very* unreachable...
            self.Spectrum.plotter._active_gui = None
            return
        for eventtype in ('button_press_event','key_press_event'):
            for key,val in iteritems(self.Spectrum.plotter.figure.canvas.callbacks.callbacks[eventtype]):
                if hasattr(val, 'func') and "event_manager" in val.func.__name__:
                    cids_to_remove.append(key)
                    if debug or self._debug: print("Removing CID #%i with attached function %s" % (key,val.func.__name__))
        for cid in cids_to_remove:
            self.Spectrum.plotter.figure.canvas.mpl_disconnect(cid)

        self.Spectrum.plotter._reconnect_matplotlib_keys()

        # Click counters - should always be reset!
        self.nclicks_b1 = 0  # button 1
        self.nclicks_b2 = 0  # button 2

        self.Spectrum.plotter._active_gui = None


    def start_interactive(self, debug=False, LoudDebug=False,
                          reset_selection=False, print_message=True,
                          clear_all_connections=True, **kwargs):
        """
        Initialize the interative session

        Parameters
        ----------
        print_message : bool
            Print the interactive help message?
        clear_all_connections : bool
            Clear all matplotlib event connections?
            (calls :func:`self.clear_all_connections`)
        reset_selection : bool
            Reset the include mask to be empty, so that you're setting up a
            fresh region.

        """
        if reset_selection:
            self.includemask[:] = False
        if print_message:
            print(self.interactive_help_message)
        if clear_all_connections:
            self.clear_all_connections()
            self.Spectrum.plotter._disconnect_matplotlib_keys()
        global_kwargs = kwargs
        def key_manager(x, *args, **kwargs):
            kwargs.update(global_kwargs)
            return self.event_manager(x, *args, debug=debug, **kwargs)
        def click_manager(x, *args, **kwargs):
            kwargs.update(global_kwargs)
            return self.event_manager(x, *args, debug=debug, **kwargs)
        key_manager.__name__ = "event_manager"
        click_manager.__name__ = "event_manager"

        self.click    = self.Spectrum.plotter.axis.figure.canvas.mpl_connect('button_press_event',click_manager)
        self.keyclick = self.Spectrum.plotter.axis.figure.canvas.mpl_connect('key_press_event',key_manager)
        self._callbacks = self.Spectrum.plotter.figure.canvas.callbacks.callbacks

        assert self._check_connections()

        self.Spectrum.plotter._active_gui = self


    def _check_connections(self, verbose=True):
        """
        Make sure the interactive session acepts user input
        """
        # check for connections
        OKclick = False
        OKkey   = False
        for cb in self._callbacks.values():
            if self.click in cb.keys():
                OKclick = True
            if self.keyclick in cb.keys():
                OKkey = True
        if self.keyclick == self.click:
            OKkey = False
        if verbose and not OKkey:
            print("Interactive session failed to connect keyboard.  Key presses will not be accepted.")
        if verbose and not OKclick:
            print("Interactive session failed to connect mouse.  Mouse clicks will not be accepted.")
        return OKkey and OKclick


    def clear_highlights(self):
        """
        Hide and remove "highlight" colors from the plot indicating the
        selected region
        """
        for p in self.button1plot:
            p.set_visible(False)
            if self.Spectrum.plotter.axis and p in self.Spectrum.plotter.axis.lines:
                if hasattr(self.Spectrum.plotter.axis.lines, 'remove'):
                    self.Spectrum.plotter.axis.lines.remove(p)
                else:
                    p.remove()
        self.button1plot = [] # I should be able to just remove from the list... but it breaks the loop...
        self.Spectrum.plotter.refresh()

    def selectregion(self, xmin=None, xmax=None, xtype='wcs', highlight=False,
                     fit_plotted_area=True, reset=False, verbose=False,
                     debug=False, use_window_limits=None, exclude=None,
                     **kwargs):
        """
        Pick a fitting region in either WCS units or pixel units

        Parameters
        ----------
        *xmin / xmax* : [ float ]
            The min/max X values to use in X-axis units (or pixel units if xtype is set).
            TAKES PRECEDENCE ALL OTHER BOOLEAN OPTIONS

        *xtype* : [ string ]
            A string specifying the xtype that xmin/xmax are specified in.  It can be either
            'wcs' or any valid xtype from :class:`pyspeckit.spectrum.units`

        *reset* : [ bool ]
            Reset the selected region to the full spectrum?  Only takes effect
            if xmin and xmax are not (both) specified.
            TAKES PRECEDENCE ALL SUBSEQUENT BOOLEAN OPTIONS

        *fit_plotted_area* : [ bool ]
            Use the plot limits *as specified in :class:`pyspeckit.spectrum.plotters`*?
            Note that this is not necessarily the same as the window plot limits!

        *use_window_limits* : [ bool ]
            Use the plot limits *as displayed*.  Defaults to self.use_window_limits
            (:attr:`pyspeckit.spectrum.interactive.use_window_limits`).
            Overwrites xmin,xmax set by plotter

        exclude: {list of length 2n,'interactive', None}
            * interactive: start an interactive session to select the
              include/exclude regions
            * list: parsed as a series of (startpoint, endpoint) in the
              spectrum's X-axis units.  Will exclude the regions between
              startpoint and endpoint
            * None: No exclusion
        """
        if debug or self._debug:
            log.info("".join(map(str, ("selectregion kwargs: ",kwargs," use_window_limits: ",use_window_limits," reset: ",reset," xmin: ",xmin, " xmax: ",xmax))))

        if reset:
            if verbose or debug or self._debug:
                print("Resetting xmin/xmax to full limits of data")
            self.xmin = 0
            # End-inclusive!
            self.xmax = self.Spectrum.data.shape[0]
            self.includemask[self.xmin:self.xmax] = True
            #raise ValueError("Need to input xmin and xmax, or have them set by plotter, for selectregion.")

        if xmin is not None and xmax is not None:
            if verbose or debug or self._debug:
                log.info("Setting xmin,xmax from keywords %g,%g" % (xmin,xmax))
            if xtype.lower() in ('wcs',) or xtype in units.xtype_dict:
                self.xmin = numpy.floor(self.Spectrum.xarr.x_to_pix(xmin))
                # End-inclusive!
                self.xmax = numpy.ceil(self.Spectrum.xarr.x_to_pix(xmax))+1
            else:
                self.xmin = xmin
                # NOT end-inclusive!  This is PYTHON indexing
                self.xmax = xmax
            self.includemask[self.xmin:self.xmax] = True
        elif (self.Spectrum.plotter.xmin is not None and
              self.Spectrum.plotter.xmax is not None and fit_plotted_area):
            if use_window_limits or (use_window_limits is None and self.use_window_limits):
                if debug or self._debug:
                    print("Resetting plotter xmin,xmax and ymin,ymax to the currently visible region")
                self.Spectrum.plotter.set_limits_from_visible_window(debug=debug)
            self.xmin = numpy.floor(self.Spectrum.xarr.x_to_pix(self.Spectrum.plotter.xmin))
            self.xmax = numpy.ceil(self.Spectrum.xarr.x_to_pix(self.Spectrum.plotter.xmax))
            if self.xmin>self.xmax:
                self.xmin,self.xmax = self.xmax,self.xmin
            # End-inclusive!  Note that this must be done after the min/max swap!
            # this feels sketchy to me, but if you don't do this the plot will not be edge-inclusive
            # that means you could do this reset operation N times to continuously shrink the plot
            self.xmax += 1
            if debug or self._debug:
                log.debug("Including all plotted area (as defined by "
                          "[plotter.xmin={0}, plotter.xmax={1}]) for "
                          "fit".format(self.Spectrum.plotter.xmin,
                                       self.Spectrum.plotter.xmax))
                log.debug("Including self.xmin:self.xmax = {0}:{1}"
                          " (and excluding the rest)".format(self.xmin,
                                                             self.xmax))
            self.includemask[self.xmin:self.xmax] = True
        else:
            if verbose:
                log.info("Left region selection unchanged."
                         "  xminpix, xmaxpix: %i,%i" % (self.xmin,self.xmax))

        if self.xmin == self.xmax:
            # Reset if there is no fitting region
            self.xmin = 0
            # End-inclusive
            self.xmax = self.Spectrum.data.shape[0]
            log.debug("Reset to full range because the endpoints were equal")
        elif self.xmin>self.xmax:
            # Swap endpoints if the axis has a negative delta-X
            self.xmin,self.xmax = self.xmax,self.xmin
            log.debug("Swapped endpoints because the left end was greater than the right")

        self.includemask[:self.xmin] = False
        self.includemask[self.xmax:] = False

        # Exclude keyword-specified excludes.  Assumes exclusion in current X array units
        log.debug("Exclude: {0}".format(exclude))
        if (isinstance(exclude, str) and (exclude == 'interactive')):
            self.start_interactive()
        elif exclude is not None and len(exclude) % 2 == 0:
            for x1,x2 in zip(exclude[::2],exclude[1::2]):
                if xtype.lower() in ('wcs',) or xtype in units.xtype_dict:
                    x1 = self.Spectrum.xarr.x_to_pix(x1)
                    # WCS units should be end-inclusive
                    x2 = self.Spectrum.xarr.x_to_pix(x2)+1
                    # correct for order if WCS units are used
                    # if pixel units are being used, we assume the user has
                    # done so intentionally
                    # TODO: if xarr, data go opposite directions, this swap
                    # doesn't work.
                    if x1 > x2:
                        x1,x2 = x2,x1
                log.debug("Exclusion pixels: {0} to {1}".format(x1,x2))
                self.includemask[x1:x2] = False
        elif exclude is not None:
            log.error("An 'exclude' keyword was specified with an odd number "
                      "of parameters, which is not permitted.")

        if highlight:
            self.highlight_fitregion()

        self._update_xminmax()

        if debug or self._debug:
            log.debug("At the end of selectregion, xmin, xmax = {0},{1}"
                      " and includemask.sum() == {2}"
                      .format(self.xmin, self.xmax, self.includemask.sum()))

    def _update_xminmax(self):
        try:
            whinclude = numpy.where(self.includemask)
            self.xmin = whinclude[0][0]
            # MUST be end-inclusive!
            self.xmax = whinclude[0][-1]+1
        except IndexError:
            pass
