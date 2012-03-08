from matplotlib.widgets import Widget,Button,Slider
from matplotlib import pyplot
import matplotlib

class FitterTool(Widget):
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
        # FIXME: The docstring seems to just abruptly end without...

        self.targetfig = targetfig
        self.specfit = specfit
        self.parlimitdict = parlimitdict

        if toolfig is None:
            tbar = matplotlib.rcParams['toolbar'] # turn off the navigation toolbar for the toolfig
            matplotlib.rcParams['toolbar'] = 'None'
            self.toolfig = pyplot.figure(figsize=(6,3))
            self.toolfig.canvas.set_window_title("Fit Sliders for "+targetfig.canvas.manager.window.title())
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
                line.set_ydata(self.specfit.get_model(line.get_xdata(),mpp))

            # update components too
            for ii,line in enumerate(self.specfit._plotted_components):
                xdata = line.get_xdata()
                modelcomponents = self.specfit.fitter.components(xdata, mpp, **self.specfit._component_kwargs)
                for jj,data in enumerate(modelcomponents):
                    if ii % 2 == jj:
                        # can have multidimensional components
                        if len(data.shape) > 1:
                            for d in (data):
                                line.set_ydata(d)
                        else:
                            line.set_ydata(data)

            self.specfit.Spectrum.plotter.refresh()


        self.sliders = []
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

            self.sliders += [Slider(ax, 
                name, vmin, vmax, valinit=value)]

            self.sliders[-1].on_changed(update)
