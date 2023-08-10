from __future__ import print_function
from six.moves import xrange
from matplotlib.widgets import Widget,Button,Slider
import matplotlib
import warnings

class dictlist(list):
    def __init__(self, *args):
        list.__init__(self, *args)

        self._dict = {}
        self._dict_index = {}
        for ii,value in enumerate(self):
            if len(value) == 2:
                self._dict[value[0]] = value[1]
                self._dict_index[value[0]] = ii
                self._dict_index[ii] = value[0]
            else:
                self._dict[ii] = value
                self._dict_index[ii] = ii

    def __getitem__(self, key):
        if type(key) is int:
            return super(dictlist,self).__getitem__(key)
        else:
            return self._dict[key]

    def __setitem__(self, key, value):
        if type(key) is int:
            super(dictlist,self).__setitem__(key,value)
            self._dict[self._dict_index[key]] = value
        else:
            if key in self._dict:
                self._dict[key] = value
                self[self._dict_index[key]] = value
            else:
                self._dict[key] = value
                self._dict_index[key] = len(self)
                self._dict_index[len(self)] = key
                self.append(value)

    def __slice__(self, s1, s2):
        pass

    def values(self):
        return [self._dict[self._dict_index[ii]] for ii in xrange(len(self))]

    def keys(self):
        return [self._dict_index[ii] for ii in xrange(len(self))]

class ModifiableSlider(Slider):

    def set_valmin(self, valmin):
        """
        Change the minimum value of the slider
        """
        self.valmin = valmin
        self.ax.set_xlim((self.valmin,self.valmax))
        if self.val < self.valmin:
            self.set_val(self.valmin)
        if self.valinit < self.valmin:
            self.valinit = (self.valmax-self.valmin)/2. + self.valmin
            if self.vline in self.ax.lines:
                if hasattr(self.ax.lines, 'remove'):
                    self.ax.lines.remove(self.vline)
                else:
                    self.vline.remove()
            self.vline = self.ax.axvline(self.valinit,0,1, color='r', lw=1)

    def set_valmax(self, valmax):
        """
        Change the maximum value of the slider
        """
        self.valmax = valmax
        self.ax.set_xlim((self.valmin,self.valmax))
        if self.val > self.valmax:
            self.set_val(self.valmax)
        if self.valinit > self.valmax:
            self.valinit = (self.valmax-self.valmin)/2. + self.valmin
            if self.vline in self.ax.lines:
                if hasattr(self.ax.lines, 'remove'):
                    self.ax.lines.remove(self.vline)
                else:
                    self.vline.remove()
            self.vline = self.ax.axvline(self.valinit,0,1, color='r', lw=1)

class FitterSliders(Widget):
    """
    A tool to adjust to subplot params of a :class:`matplotlib.figure.Figure`
    """
    def __init__(self, specfit, targetfig, npars=1, toolfig=None, parlimitdict={}):
        """
        *targetfig*
            The figure instance to adjust

        *toolfig*
            The figure instance to embed the subplot tool into. If
            None, a default figure will be created. If you are using
            this from the GUI
        """

        self.targetfig = targetfig
        self.specfit = specfit
        self.parlimitdict = parlimitdict

        from matplotlib import pyplot

        if toolfig is None:
            tbar = matplotlib.rcParams['toolbar'] # turn off the navigation toolbar for the toolfig
            matplotlib.rcParams['toolbar'] = 'None'
            self.toolfig = pyplot.figure(figsize=(6,3))
            if hasattr(targetfig.canvas.manager,'window'):
                if hasattr(targetfig.canvas.manager.window, 'title') and hasattr(targetfig.canvas, 'set_window_title'):
                    self.toolfig.canvas.set_window_title("Fit Sliders for "+targetfig.canvas.manager.window.title())
                elif hasattr(targetfig.canvas.manager.window, 'windowTitle') and hasattr(targetfig.canvas, 'set_window_title'):
                    self.toolfig.canvas.set_window_title("Fit Sliders for "+targetfig.canvas.manager.window.windowTitle())
                else:
                    warnings.warn("Only Qt4 and TkAgg support window titles (apparently)")
            self.toolfig.subplots_adjust(top=0.9,left=0.2,right=0.9)
            matplotlib.rcParams['toolbar'] = tbar
        else:
            self.toolfig = toolfig
            self.toolfig.subplots_adjust(left=0.2, right=0.9)


        bax = self.toolfig.add_axes([0.8, 0.05, 0.15, 0.075])
        self.buttonreset = Button(bax, 'Reset')



        self.set_sliders(parlimitdict)

        def reset(event):
            thisdrawon = self.drawon

            self.drawon = False

            # store the drawon state of each slider
            bs = []
            for slider in self.sliders:
                bs.append(slider.drawon)
                slider.drawon = False

            # reset the slider to the initial position
            for slider in self.sliders:
                slider.reset()

            # reset drawon
            for slider, b in zip(self.sliders, bs):
                slider.drawon = b

            # draw the canvas
            self.drawon = thisdrawon
            if self.drawon:
                self.toolfig.canvas.draw()
                self.targetfig.canvas.draw()


        # during reset there can be a temporary invalid state
        # depending on the order of the reset so we turn off
        # validation for the resetting
        if hasattr(self.toolfig.subplotpars, 'validate'):
            validate = self.toolfig.subplotpars.validate
            self.toolfig.subplotpars.validate = False
            self.buttonreset.on_clicked(reset)
            self.toolfig.subplotpars.validate = validate


    def clear_sliders(self):
        """
        Get rid of the sliders...
        """
        try:
            for sl in self.sliders:
                sl.ax.remove()
        except NotImplementedError:
            for sl in self.sliders:
                self.specfit.Spectrum.plotter.figure.delaxes(sl.ax)

        self.specfit.Spectrum.plotter.refresh()

    def set_sliders(self, parlimitdict={}):
        """
        Set the slider properties, actions, and values

        can also reset their limits
        """

        def update(value):
            mpp = [slider.val for slider in self.sliders]
            for line in self.specfit.modelplot:
                line.set_ydata(self.specfit.get_model_frompars(line.get_xdata(),mpp))

            # update components too
            for ii,line in enumerate(self.specfit._plotted_components):
                xdata = line.get_xdata()
                modelcomponents = self.specfit.fitter.components(xdata,
                                                                 mpp,
                                                                 **self.specfit._component_kwargs)
                for jj,data in enumerate(modelcomponents):
                    if ii % 2 == jj:
                        # can have multidimensional components
                        if len(data.shape) > 1:
                            for d in (data):
                                line.set_ydata(d)
                        else:
                            line.set_ydata(data)

            self.specfit.Spectrum.plotter.refresh()


        self.sliders = dictlist()
        npars = len(self.specfit.parinfo)
        for param in self.specfit.parinfo:
            name = param['parname']
            value = param['value']
            limited = param['limited']
            limits = param['limits']

            # make one less subplot so that there's room for buttons
            # param['n'] is zero-indexed, subplots are 1-indexed
            ax = self.toolfig.add_subplot(npars+1,1,param['n']+1)
            ax.set_navigate(False)

            if name in parlimitdict:
                limits = parlimitdict[name]
                limited = [True,True]
            if limited[0]:
                vmin = limits[0]
            elif value != 0:
                vmin = min([value/4.0,value*4.0])
            else:
                vmin = -1

            if limited[1]:
                vmax = limits[1]
            elif value != 0:
                vmax = max([value/4.0,value*4.0])
            else:
                vmax = 1
            try:
                self.sliders[name] = ModifiableSlider(ax,
                    name, vmin, vmax, valinit=value)
            except ValueError:
                self.sliders[name] = ModifiableSlider(ax,
                    name, vmin.value, vmax.value, valinit=value)

            self.sliders[-1].on_changed(update)



    def get_values(self):
        return [s.val for s in self.sliders]

class FitterTools(Widget):
    """
    A tool to monitor and play with :class:`pyspeckit.spectrum.fitter` properties

--------------------------
| Baseline range  [x,x]  |
| Baseline order  -      |
| (Baseline subtracted)  |
|                        |
| Fitter range    [x,x]  |
| Fitter type    ------- |
| Fitter Guesses  [p,w]  |
|        ...      ...    |
|        ...      ...    |
|                        |
| (Fit) (BL fit) (reset) |
--------------------------


    """
    def __init__(self, specfit, targetfig, toolfig=None, nsubplots=12):
        """
        *targetfig*
            The figure instance to adjust

        *toolfig*
            The figure instance to embed the subplot tool into. If
            None, a default figure will be created. If you are using
            this from the GUI
        """

        self.targetfig = targetfig
        self.specfit = specfit
        self.baseline = specfit.Spectrum.baseline
        self.plotter = specfit.Spectrum.plotter

        from matplotlib import pyplot

        if toolfig is None:
            tbar = matplotlib.rcParams['toolbar'] # turn off the navigation toolbar for the toolfig
            matplotlib.rcParams['toolbar'] = 'None'
            self.toolfig = pyplot.figure(figsize=(6,3))
            if hasattr(self.toolfig.canvas, 'set_window_title'):
                self.toolfig.canvas.set_window_title("Fit Tools for "+targetfig.canvas.manager.window.title())
            self.toolfig.subplots_adjust(top=0.9,left=0.05,right=0.95)
            matplotlib.rcParams['toolbar'] = tbar
        else:
            self.toolfig = toolfig
            self.toolfig.subplots_adjust(left=0.0, right=1.0)


        #bax = self.toolfig.add_axes([0.6, 0.05, 0.15, 0.075])
        #self.buttonrefresh = Button(bax, 'Refresh')

        # buttons ruin everything.
        # fax = self.toolfig.add_axes([0.1, 0.05, 0.15, 0.075])
        # self.buttonfit = Button(fax, 'Fit')
        #
        # resetax = self.toolfig.add_axes([0.7, 0.05, 0.15, 0.075])
        # self.buttonreset = Button(resetax, 'Reset')

        # resetblax = self.toolfig.add_axes([0.3, 0.05, 0.15, 0.075])
        # self.buttonresetbl = Button(resetblax, 'Reset BL')

        # resetfitax = self.toolfig.add_axes([0.5, 0.05, 0.15, 0.075])
        # self.buttonresetfit = Button(resetfitax, 'Reset fit')

        def refresh(event):
            thisdrawon = self.drawon

            self.drawon = False

            self.update_information()

            # draw the canvas
            self.drawon = thisdrawon
            if self.drawon:
                self.toolfig.canvas.draw()
                self.targetfig.canvas.draw()

        def fit(event):
            self.specfit.button3action(event)

        def reset_fit(event):
            self.specfit.guesses = []
            self.specfit.npeaks = 0
            self.specfit.includemask[:] = True
            self.refresh(event)

        def reset_baseline(event):
            self.baseline.unsubtract()
            self.refresh(event)

        def reset(event):
            reset_baseline(event)
            reset_fit(event)
            self.plotter()
            self.refresh(event)

        # during refresh there can be a temporary invalid state
        # depending on the order of the refresh so we turn off
        # validation for the refreshting
        #validate = self.toolfig.subplotpars.validate
        #self.toolfig.subplotpars.validate = False
        #self.buttonrefresh.on_clicked(refresh)
        #self.toolfig.subplotpars.validate = validate

        # these break everything.
        # self.buttonfit.on_clicked(fit)
        # self.buttonresetfit.on_clicked(reset_fit)
        # self.buttonresetbl.on_clicked(reset_baseline)
        # self.buttonreset.on_clicked(reset)


        #menuitems = []
        #for label in ('polynomial','blackbody','log-poly'):
        #    def on_select(item):
        #        print 'you selected', item.labelstr
        #    item = MenuItem(fig, label, props=props, hoverprops=hoverprops,
        #                    on_select=on_select)
        #    menuitems.append(item)

        #menu = Menu(fig, menuitems)


        self.axes = [self.toolfig.add_subplot(nsubplots,1,spnum, frame_on=False, navigate=False, xticks=[], yticks=[])
                for spnum in xrange(1,nsubplots+1)]
        #self.axes = self.toolfig.add_axes([0,0,1,1])

        self.use_axes = [0,1,2,4,5,6,7,8,9,10,11]
        self.labels = dict([(axnum,None) for axnum in self.use_axes])
        self.update_information()

        self.targetfig.canvas.mpl_connect('button_press_event',self.refresh)
        self.targetfig.canvas.mpl_connect('key_press_event',self.refresh)
        self.targetfig.canvas.mpl_connect('draw_event',self.refresh)

    def refresh(self, event):
        try:
            thisdrawon = self.drawon

            self.drawon = False

            self.update_information()

            # draw the canvas
            self.drawon = thisdrawon
            if self.drawon:
                self.toolfig.canvas.draw()
        except:
            # ALWAYS fail silently
            # this is TERRIBLE coding practice, but I have no idea how to tell the object to disconnect
            # when the figure is closed
            pass

    def update_information(self, **kwargs):
        self.information = [
            ("Baseline Range","(%g,%g)" % (self.baseline.xmin,self.baseline.xmax)),
            ("Baseline Order","%i" % (self.baseline.order)),
            ("Baseline Subtracted?","%s" % (self.baseline.subtracted)),
            ("Fitter Range","(%g,%g)" % (self.specfit.xmin,self.specfit.xmax)),
            ("Fitter Type","%s" % (self.specfit.fittype)),
            ]

        for ii in xrange(self.specfit.npeaks):
            guesses = tuple(self.specfit.guesses[ii:ii+3])
            if len(guesses) == 3:
                self.information += [("Fitter guesses%i:" % ii , "p: %g c: %g w: %g" % guesses) ]
            else:
                break

        self.show_labels(**kwargs)

    def show_selected_region(self):
        self.specfit.highlight_fitregion()

    def show_label(self, axis, text, xloc=0.0, yloc=0.5, **kwargs):
        return axis.text(xloc, yloc, text, **kwargs)

    def show_value(self, axis, text, xloc=0.5, yloc=0.5, **kwargs):
        return axis.text(xloc, yloc, text, **kwargs)

    def show_labels(self, **kwargs):
        for axnum,(label,text) in zip(self.use_axes, self.information):
            if self.labels[axnum] is not None and len(self.labels[axnum]) == 2:
                labelobject,textobject = self.labels[axnum]
                labelobject.set_label(label)
                textobject.set_text(text)
            else:
                self.labels[axnum] = (self.show_label(self.axes[axnum],label),
                        self.show_value(self.axes[axnum],text))

    def update_info_texts(self):
        for newtext,textobject in zip(self.information.values(), self.info_texts):
            textobject.set_text(newtext)


#import parinfo
#
#class ParameterButton(parinfo.Parinfo):
#    """
#    A class to manipulate individual parameter values
#    """
#    def __init__(self,
