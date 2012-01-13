"""
==================
Window Interaction
==================

A general module for selecting regions and inputting guesses via the
interactive window.
"""
import numpy

class Interactive(object):

    def __init__(Spectrum, guesses=None,
            interactive_help_message="Replace this message"):
        """
        Declare interactive variables.

        Must have a parent Spectrum class

        **Must declare button2action and button3action**
        """
        self.Spectrum = Spectrum

        self.interactive_help_message = interactive_help_message
        self.includemask = None
        self.excludemask = None
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

    def event_manager(self, event):
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

            if debug:
                print "button: ",button

            if button in ('p','P','1',1,'i','a'):
                # button one is always region selection
                self.selectregion_interactive(event,debug=debug)
            elif button in ('e','x','r','E','X','R'):
                # exclude/delete/remove
                self.selectregion_interactive(event, mark_include=False, debug=debug)
            elif button in ('m','M','2',2):
                self.button2action(event,debug=debug)
            elif button in ('d','D','3',3):
                self.button3action(event,debug=debug)
            elif button in ('?'):
                print self.interactive_help_message
            elif button in self.Registry.fitkeys:
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

            self.button1plot += self.Spectrum.plotter.axis.plot(
                    self.Spectrum.xarr[self._xclick1:self._xclick2],
                    self.Spectrum.data[self._xclick1:self._xclick2]+self.Spectrum.plotter.offset,
                    drawstyle='steps-mid',
                    color='g',alpha=0.5)
            self.Spectrum.plotter.refresh()
            self.includemask[self._xclick1:self._xclick2] = mark_include
            if debug: print "Click 2: clickx=%i xmin=%i, xmax=%i" % (xpix,self.xmin,self.xmax)

    def show_fitregion(self):
        bad = self.Spectrum.data*0
        bad[self.excludemask] = np.nan
        self.button1plot += self.Spectrum.plotter.axis.plot(
                self.Spectrum.xarr,
                self.Spectrum.data+self.Spectrum.plotter.offset+bad,
                drawstyle='steps-mid',
                color='g',alpha=0.5)
        self.Spectrum.plotter.refresh()

    def firstclick_selection(self, include_all=False):
        """
        Initialize the include/exclude mask
        """

        if include_all:
            # default to including everything
            self.includemask = numpy.array(self.Spectrum.data, dtype='bool') + True
        else:
            # default to including nothing
            self.includemask = numpy.array(self.Spectrum.data, dtype='bool') * False

    def eventhandler_debug(self,event):
        return self.makeguess(event,debug=True)

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
            self.button2plot += [self.specplotter.axis.scatter(event.xdata,event.ydata,marker='x',c='r')]
            #self.specplotter.refresh() #plot(**self.specplotter.plotkwargs)
        elif self.nclicks_b2 % 2 == 1:
            self.guesses[-1] = abs(event.xdata-self.guesses[-2]) / np.sqrt(2*np.log(2))
            self.nclicks_b2 += 1
            if debug: print "Width %i click %i at x,y %g,%g" % (self.npeaks,self.nclicks_b2,event.xdata,event.ydata)
            self.button2plot += self.specplotter.axis.plot([event.xdata,
                2*self.guesses[-2]-event.xdata],[event.ydata]*2,
                color='r')
            #self.specplotter.refresh() #plot(**self.specplotter.plotkwargs)
            if self.auto:
                self.auto = False
            if self.nclicks_b2 / 2 > self.npeaks:
                print "There have been %i middle-clicks but there are only %i gaussians" % (self.nclicks_b2,self.npeaks)
                self.npeaks += 1

    def firstclick_guess(self):
        """
        Initialize self.guesses
        """
        if self.guesses is None:
            self.guesses = []
        elif len(self.guesses) > 0:
            for ii in xrange(len(guesses)): 
                self.guesses.pop()

    def button2_baseline(self, event, debug):
        self.Spectrum.plotter.figure.canvas.mpl_disconnect(self.click)
        self.dofit(include=self.includepix, includeunits='pix', **kwargs)
        for p in self.button1plot
            p.set_visible(False)
            if p in self.Spectrum.plotter.axis.lines: self.Spectrum.plotter.axis.lines.remove(p)
        self.button1plot=[] # I should be able to just remove from the list... but it breaks the loop...
        self.Spectrum.plotter.refresh()
        if debug: print "Click to fit.  Includepix: %s" % (str(self.includepix))

    def clear_all_connections(self):
        """
        Prevent overlapping interactive sessions
        """
        # this is really ugly, but needs to be done in order to prevent multiple overlapping calls...
        cids_to_remove = []
        for eventtype in ('button_press_event','key_press_event'):
            for key,val in self.Spectrum.plotter.figure.canvas.callbacks.callbacks[eventtype].iteritems():
                if "event_handler" in val.func.__name__:
                    cids_to_remove.append(key)
                    if debug: print "Removing CID #%i with attached function %s" % (key,val.func.__name__)
        for cid in cids_to_remove:
            self.Spectrum.plotter.figure.canvas.mpl_disconnect(cid)


