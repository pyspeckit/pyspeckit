"""
==================
Window Interaction
==================

A general module for selecting regions and inputting guesses via the
interactive window.
"""
import numpy

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
        self.includemask = numpy.array( self.Spectrum.data.astype('bool') + True )
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

    def event_manager(self, event, debug=False):
        """
        Decide what to do given input (click, keypress, etc.)
        """
        toolbar = self.Spectrum.plotter.figure.canvas.manager.toolbar
        self.event_history.append(event)

        if toolbar.mode == '' and self.Spectrum.plotter.axis in event.canvas.figure.axes:
            if hasattr(event,'button'):
                button = event.button
            elif hasattr(event,'key'):
                button = event.key

            if event.xdata is None or event.ydata is None:
                return

            if debug:
                print "button: ",button," x,y: ",event.xdata,event.ydata

            if button in ('p','P','1',1,'i','a'): # p for... parea?  a for area.  i for include
                # button one is always region selection
                self.selectregion_interactive(event,debug=debug)
            elif button in ('e','x','r','E','X','R'): # e for exclude, x for x-clude, r for remove
                # exclude/delete/remove
                self.selectregion_interactive(event, mark_include=False, debug=debug)
            elif button in ('m','M','2',2): # m for mark
                if debug: print "Button 2 action"
                self.button2action(event,debug=debug)
            elif button in ('d','D','3',3): # d for done
                self.button3action(event,debug=debug)
            elif button in ('?'):
                print self.interactive_help_message
            elif hasattr(self,'Registry') and button in self.Registry.fitkeys:
                fittername = self.Registry.fitkeys[button]
                if fittername in self.Registry.multifitters:
                    self.fitter = self.Registry.multifitters[fittername]
                    self.fittype = fittername
                    print "Selected multi-fitter %s" % fittername
                elif fittername in self.Registry.singlefitters:
                    self.fitter = self.Registry.singlefitters[fittername]
                    self.fittype = fittername
                    print "Selected single-fitter %s" % fittername
                else: 
                    print "ERROR: Did not find fitter %s" % fittername
            if self.Spectrum.plotter.autorefresh: self.Spectrum.plotter.refresh()


    def selectregion_interactive(self, event, mark_include=True, debug=False, **kwargs):
        """
        select regions for baseline fitting
        """

        xpix = self.Spectrum.xarr.x_to_pix(event.xdata)
        if self.xclicks == []:
            self.firstclick_selection(not mark_include)

        if self.nclicks_b1 == 0:
            self.nclicks_b1 = 1
            self._xclick1 = xpix
            self.xclicks.append(xpix)
            if debug: print "Click 1: clickx=%i xmin=%i, xmax=%i" % (xpix,self.xmin,self.xmax)
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
        
            if mark_include: 
                color = 'g'
                linewidth = 0.5
            else:
                try:
                    color = self.Spectrum.plotter._spectrumplot[0].get_color()
                    linewidth = self.Spectrum.plotter._spectrumplot[0].get_linewidth()
                except IndexError:
                    linewidth = 0.5
                    color = 'k'

            self.button1plot += self.Spectrum.plotter.axis.plot(
                    self.Spectrum.xarr[self._xclick1:self._xclick2],
                    self.Spectrum.data[self._xclick1:self._xclick2]+self.Spectrum.plotter.offset,
                    drawstyle='steps-mid',
                    color=color,
                    linewidth=linewidth)
            self.Spectrum.plotter.refresh()
            self.includemask[self._xclick1:self._xclick2] = mark_include
            if debug: print "Click 2: clickx=%i xmin=%i, xmax=%i" % (xpix,self.xmin,self.xmax)

    def highlight_fitregion(self,  drawstyle='steps-mid', color='g',
            clear_highlights=True, **kwargs):
        """
        Re-highlight the fitted region

        kwargs are passed to `matplotlib.plot`
        """
        if clear_highlights: self.clear_highlights()
        bad = self.Spectrum.data*0
        bad[True-self.includemask] = numpy.nan
        self.button1plot += self.Spectrum.plotter.axis.plot(
                self.Spectrum.xarr,
                self.Spectrum.data+self.Spectrum.plotter.offset+bad,
                drawstyle=drawstyle, color=color, 
                **kwargs)
        self.Spectrum.plotter.refresh()

    def firstclick_selection(self, include_all=False):
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

    def guesspeakwidth(self,event,debug=False):
        """
        Interactively guess the peak height and width from user input

        Width is assumed to be half-width-half-max
        """
        if self.nclicks_b2 == 0:
            self.firstclick_guess()
        if self.nclicks_b2 % 2 == 0:
            self.guesses += [event.ydata,event.xdata,1]
            self.npeaks += 1
            self.nclicks_b2 += 1
            if debug: print "Peak %i click %i at x,y %g,%g" % (self.npeaks,self.nclicks_b2,event.xdata,event.ydata)
            self.button2plot += [self.Spectrum.plotter.axis.scatter(event.xdata,event.ydata,marker='x',c='r')]
            #self.Spectrum.plotter.refresh() #plot(**self.Spectrum.plotter.plotkwargs)
        elif self.nclicks_b2 % 2 == 1:
            self.guesses[-1] = abs(event.xdata-self.guesses[-2]) / numpy.sqrt(2*numpy.log(2))
            self.nclicks_b2 += 1
            if debug: print "Width %i click %i at x,y %g,%g" % (self.npeaks,self.nclicks_b2,event.xdata,event.ydata)
            self.button2plot += self.Spectrum.plotter.axis.plot([event.xdata,
                2*self.guesses[-2]-event.xdata],[event.ydata]*2,
                color='r')
            #self.Spectrum.plotter.refresh() #plot(**self.Spectrum.plotter.plotkwargs)
            if self.nclicks_b2 / 2 > self.npeaks:
                print "There have been %i middle-clicks but there are only %i gaussians" % (self.nclicks_b2,self.npeaks)
                self.npeaks += 1

    def firstclick_guess(self):
        """
        Initialize self.guesses
        """
        self.Spectrum.plotter.axis.set_autoscale_on(False)
        if self.guesses is None:
            self.guesses = []
        elif len(self.guesses) > 0:
            for ii in xrange(len(guesses)): 
                self.guesses.pop()

    def clear_all_connections(self, debug=False):
        """
        Prevent overlapping interactive sessions
        """
        # this is really ugly, but needs to be done in order to prevent multiple overlapping calls...
        cids_to_remove = []
        if not hasattr(self.Spectrum.plotter.figure,'canvas'):
            # just quit out; saves a tab...
            if debug: print "Didn't find a canvas, quitting."
            return
        for eventtype in ('button_press_event','key_press_event'):
            for key,val in self.Spectrum.plotter.figure.canvas.callbacks.callbacks[eventtype].iteritems():
                if "event_manager" in val.func.__name__:
                    cids_to_remove.append(key)
                    if debug: print "Removing CID #%i with attached function %s" % (key,val.func.__name__)
        for cid in cids_to_remove:
            self.Spectrum.plotter.figure.canvas.mpl_disconnect(cid)


    def start_interactive(self, debug=False, LoudDebug=False, print_message=True, clear_all_connections=True, **kwargs):
        """

        """
        if print_message: print self.interactive_help_message
        if clear_all_connections: self.clear_all_connections()
        event_manager = lambda(x): self.event_manager(x, debug=debug, **kwargs)
        event_manager.__name__ = "event_manager"
        self.click = self.Spectrum.plotter.axis.figure.canvas.mpl_connect('button_press_event',event_manager)
        self.keyclick = self.Spectrum.plotter.axis.figure.canvas.mpl_connect('key_press_event',event_manager)
        self._callbacks = self.Spectrum.plotter.figure.canvas.callbacks.callbacks

    def clear_highlights(self):
        for p in self.button1plot:
            p.set_visible(False)
            if p in self.Spectrum.plotter.axis.lines: self.Spectrum.plotter.axis.lines.remove(p)
        self.button1plot=[] # I should be able to just remove from the list... but it breaks the loop...
        self.Spectrum.plotter.refresh()

    def selectregion(self, xmin=None, xmax=None, xtype='wcs', highlight=False,
            fit_plotted_area=True, reset=False, verbose=False, debug=False,
            **kwargs):
        """
        Pick a fitting region in either WCS units or pixel units
        """
        if xmin is not None and xmax is not None:
            if xtype in ('wcs','WCS','velo','velocity','wavelength','frequency','freq','wav'):
                self.xmin = self.Spectrum.xarr.x_to_pix(xmin)
                self.xmax = self.Spectrum.xarr.x_to_pix(xmax)
            else:
                self.xmin = xmin
                self.xmax = xmax
        elif self.Spectrum.plotter.xmin is not None and self.Spectrum.plotter.xmax is not None and fit_plotted_area:
            self.xmin = self.Spectrum.xarr.x_to_pix(self.Spectrum.plotter.xmin)
            self.xmax = self.Spectrum.xarr.x_to_pix(self.Spectrum.plotter.xmax)
        elif reset:
            self.xmin = 0
            self.xmax = self.Spectrum.data.shape[0]
            #raise ValueError("Need to input xmin and xmax, or have them set by plotter, for selectregion.")
        else:
            if verbose: print "Left region selection unchanged.  xminpix, xmaxpix: %i,%i" % (self.xmin,self.xmax)
        
        if self.xmin == self.xmax:
            # Reset if there is no fitting region
            self.xmin = 0
            self.xmax = self.Spectrum.data.shape[0]
            if debug: print "Reset to full range because the endpoints were equal"
        elif self.xmin>self.xmax: 
            # Swap endpoints if the axis has a negative delta-X
            self.xmin,self.xmax = self.xmax,self.xmin
            if debug: print "Swapped endpoints because the left end was greater than the right"

        self.includemask[:self.xmin] = False
        self.includemask[self.xmax:] = False

        if highlight:
            self.highlight_fitregion()
