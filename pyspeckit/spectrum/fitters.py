from __future__ import print_function
import matplotlib
import numpy as np
import copy
import re
import warnings

from astropy import log
from astropy import units as u
from six.moves import xrange
from six import string_types

from ..config import mycfg
from ..config import ConfigDescriptor as cfgdec
from . import units
from . import models
from ..specwarnings import warn
from . import interactive
from . import history
from . import widgets

from pyspeckit.spectrum.units import SpectroscopicAxis

class Registry(object):
    """
    This class is a simple wrapper to prevent fitter properties from being globals
    """

    def __init__(self):
        self.npars = {}
        self.multifitters = {}
        #to delete
        self.peakbgfitters = {}
        self.fitkeys = {}
        self.associatedkeys = {}

        self._interactive_help_message_root = """

'?' will print this help message again. The keys denoted by surrounding / / are
mnemonics.
1. Left-click or hit 'p' (/p/ick) with the cursor over the plot at both of the
two desired X-values to select a fitting range.  You can e/x/clude parts of the
spectrum by hitting 'x' at two positions.
2. Then /m/iddle-click or hit 'm' twice to select (/m/ark) a peak and width -
the first mark should be on the peak of the line, the second should be at the
approximate half-max point on the curve.
3. When you're done, right-click or hit 'd' to perform the fit and disconnect
the mouse and keyboard (/d/isconnect because you're /d/one).  Any time before
you're /d/one, you can select a different fitter (see below).

To /c/ancel or /c/lear all connections, press 'c'

'?' : get help (this message)
'c' : cancel / clear
'p','1' : pick / selection region for fitting
'm','2' : mark / identify a peak
'd','3' : done / do the fit, then disconnect the fitter
'i' : individual components / show each fitted component

You can select different fitters to use with the interactive fitting routine.
The default is gaussian ('g'), all options are listed below:
        """
        self._make_interactive_help_message()

    def __copy__(self):
        # http://stackoverflow.com/questions/1500718/what-is-the-right-way-to-override-the-copy-deepcopy-operations-on-an-object-in-p
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, copy.deepcopy(v, memo))
        return result

    def add_fitter(self, name, function, npars, override=False, key=None,
                   multisingle=None):
        '''
        Register a fitter function.

        Parameters
        ----------
        name: string
            The fit function name.
        function: function
            The fitter function.  Single-fitters should take npars + 1 input
            parameters, where the +1 is for a 0th order baseline fit.  They
            should accept an X-axis and data and standard fitting-function
            inputs (see, e.g., gaussfitter).  Multi-fitters should take N *
            npars, but should also operate on X-axis and data arguments.
        npars: int
            How many parameters does the function being fit accept?

        Other Parameters
        ----------------
        override: True | False
            Whether to override any existing type if already present.
        key: char
            Key to select the fitter in interactive mode
        '''
        if multisingle is not None:
            warn("The 'multisingle' keyword is no longer required.",
                 DeprecationWarning)

        if not name in self.peakbgfitters or override:
            self.peakbgfitters[name] = function

        if not name in self.multifitters or override:
            self.multifitters[name] = function

        if key is not None:
            self.fitkeys[key] = name
            self._make_interactive_help_message()
        self.npars[name] = npars
        self.associated_keys = dict(zip(self.fitkeys.values(),self.fitkeys.keys()))

    def _make_interactive_help_message(self):
        """
        Generate the interactive help message from the fitkeys
        """
        self.interactive_help_message = (
                self._interactive_help_message_root +
                "\n" +
                "\n".join(["'%s' - select fitter %s" % (key,name) for key,name in self.fitkeys.items()]) +
                "\n" # trailing \n so that users' input is on a fresh line
                )


# Declare default registry built in for all spectra
default_Registry = Registry()
default_Registry.add_fitter('ammonia',models.ammonia_model(),6,key='a')
default_Registry.add_fitter('cold_ammonia',models.ammonia.cold_ammonia_model(),6)
default_Registry.add_fitter('ammonia_tau',models.ammonia_model_vtau(),6)
# not implemented default_Registry.add_fitter(Registry,'ammonia',models.ammonia_model( ),6, ,key='A')
default_Registry.add_fitter('formaldehyde',models.formaldehyde_fitter,3,key='F') # CAN'T USE f!  reserved for fitting
# do'nt override default_Registry.add_fitter('formaldehyde',models.formaldehyde_vheight_fitter,3)
default_Registry.add_fitter('gaussian',models.gaussian_fitter(),3,key='g')
default_Registry.add_fitter('vheightgaussian',models.gaussian_vheight_fitter(),4)
default_Registry.add_fitter('voigt',models.voigt_fitter(),4,key='v')
default_Registry.add_fitter('lorentzian',models.lorentzian_fitter(),3,key='L')
#default_Registry.add_fitter('hill5',models.hill5infall.hill5_fitter,5)
#default_Registry.add_fitter('hcn',models.hcn.hcn_vtau_fitter,4)


class Specfit(interactive.Interactive):

    def __init__(self, Spectrum, Registry=None):
        super(Specfit, self).__init__(Spectrum,
                                      interactive_help_message=Registry.interactive_help_message)
        self.model = None
        self.parinfo = None
        self.modelpars = None
        self.modelerrs = None
        self.modelplot = []
        self.modelcomponents = None
        self._plotted_components = []
        self.npeaks = 0
        #self.nclicks_b1 = 0
        #self.nclicks_b2 = 0
        #self.xmin = 0
        #self.xmax = Spectrum.data.shape[0]
        self.button2action = self.guesspeakwidth
        self.guesses = []
        self.click = 0
        self.fitkwargs = {}
        self.auto = False
        self.fitleg=None
        self.residuals=None
        self.setfitspec()
        self.fittype = 'gaussian'
        self.measurements = None
        self.vheight=False # vheight must be a boolean, can't be none
        self._component_kwargs = {}
        self.Registry = Registry
        self.autoannotate = mycfg['autoannotate']
        self.EQW_plots = []
        #self.seterrspec()

    @property
    def fitter(self):
        if hasattr(self, '_fitter'):
            return self._fitter
        else:
            raise AttributeError("The 'specfit' object has no 'fitter' yet.  "
                                 "This means you haven't yet run a fit.  The "
                                 "fitter is not accessible until after a fit "
                                 "has been run.")

    @fitter.setter
    def fitter(self, value):
        self._fitter = value

    @cfgdec
    def __call__(self, interactive=False, usemoments=True,
                 fittype=None,
                 clear_all_connections=True, debug=False, guesses='moments',
                 parinfo=None, save=True, annotate=None, show_components=None,
                 use_lmfit=False, verbose=True, clear=True,
                 reset_selection=True,
                 fit_plotted_area=True, use_window_limits=None, vheight=None,
                 exclude=None, **kwargs):
        """
        Fit model functions to a spectrum

        Parameters
        ----------
        interactive : boolean
            The plotter window will go into interactive mode.  See
            self.interactive_help_message for details on how to use the
            interactive fitter.
        fittype : str
            [passed to fitting codes; defaults to gaussian]
            The model to use.  Model must be registered in self.Registry.
            gaussian, lorentzian, and voigt profiles are registered by default
        guesses : list or 'moments'
            A list of guesses.  Guesses must have length = n*number of parameters
            in model.  Guesses are *required* for multifit fits (there is no
            automated guessing for most models)
            EXAMPLE: for single-fit gaussian
            guesses = [height,amplitude,center,width]
            for multi-fit gaussian, it is
            [amplitude, center, width]
            You can also pass the keyword string 'moments' to have the moments
            be used to automatically determine the guesses for a *single* peak
        parinfo : `pyspeckit.spectrum.parinfo.ParinfoList`
            An alternative way to specify guesses.  Supercedes guesses.
        use_lmfit : boolean
            If lmfit-py (https://github.com/newville/lmfit-py) is installed, you
            can use it instead of the pure-python (but slow) mpfit.
        reset_selection : boolean
            Override any selections previously made using `fit_plotted_area` or
            other keywords?
        fit_plotted_area : boolean
            If no other limits are specified, the plotter's xmin/xmax will be
            used to define the fit region.  Only respects the x-axis limits,
            not the y-axis limits.
        use_window_limits : boolean
            If ``fit_plotted_area==True`` and no other limits are specified,
            will use the displayed window area (as set by the zoom tools) as
            the fitting range.  Only respects the x-axis limits, not the y-axis
            limits.
        exclude : None or list
            Passed to selectregion; specifies regions to exclude in xarr units


        Plotter-related Parameters
        --------------------------
        annotate : None or boolean
            If None, will use the default stored in self.annotate, otherwise
            overwrites.  Annotations will appear on the plot if a plotter
            exists.
        show_components : boolean
            Show the individual components of a multi-component fit (defaults
            to blue)
        clear : boolean
            Clear previous fitter plots before overplotting the fit?


        Advanced Parameters
        -------------------
        clear_all_connections : boolean
            Clear all of the interactive connections from a previous interactive
            session (e.g., a baseline fitting session) before continuing?
        usemoments : boolean
            Use the moments of the spectrum as input guesses.  Only works
            for gaussian and gaussian-like models.  Only works for single-fit
            mode (not multifit)
            DEPRECATED
        debug : boolean
            Print debug statements?
        save : boolean
            Save best fits in the FITS header as keywords?  ONLY IMPLEMENTED
            FOR GAUSSIANS
        verbose : boolean
            Print out extra stuff
        vheight : None or boolean
            if None, defaults to self.vheight, otherwise overrides
            Determines whether a 0th order baseline will be fit along with the
            line

        """
        if clear:
            self.clear()
        if reset_selection:
            self.selectregion(verbose=verbose, debug=debug,
                              fit_plotted_area=fit_plotted_area,
                              exclude=exclude,
                              use_window_limits=use_window_limits, **kwargs)
        for arg in ['xmin','xmax','xtype','reset']:
            if arg in kwargs:
                kwargs.pop(arg)

        if fittype is not None:
            self.fittype = fittype
            self.fitter = self.Registry.multifitters[self.fittype]

        if 'multifit' in kwargs:
            kwargs.pop('multifit')
            log.warning("The multifit keyword is no longer required.  All fits "
                        "allow for multiple components.", DeprecationWarning)

        if 'guess' in kwargs:
            if guesses is None:
                guesses = kwargs.pop('guess')
                log.warning("The keyword 'guess' is nonstandard; please use 'guesses'")
            else:
                raise ValueError("Received keywords 'guess' and 'guesses'.  "
                                 "Please only use 'guesses'")

        self.npeaks = 0
        self.fitkwargs = kwargs
        log.debug("Additional keyword arguments passed to fitter are: {0}"
                  .format(kwargs))
        if interactive:
            if self.Spectrum.plotter.axis is None:
                raise Exception("Interactive fitting requires a plotter.")
            # reset button count & guesses on every __call__
            self.nclicks_b1 = 0
            self.nclicks_b2 = 0
            self.guesses = []

            self.start_interactive(clear_all_connections=clear_all_connections,
                                   reset_selection=True,
                                   debug=debug, **kwargs)
        elif (self.fittype in self.Registry.multifitters
              or guesses is not None
              or parinfo is not None):
            if guesses is None and parinfo is None:
                raise ValueError("You must input guesses when using multifit."
                                 "  Also, baseline (continuum fit) first!")
            elif parinfo is not None:
                self.guesses = parinfo.values
                self.parinfo = parinfo
                self.multifit(show_components=show_components, verbose=verbose,
                              debug=debug, use_lmfit=use_lmfit,
                              annotate=annotate, parinfo=parinfo,
                              guesses=None, **kwargs)
            elif guesses is not None:
                if isinstance(guesses, tuple):
                    guesses = list(guesses)
                self.guesses = guesses
                self.multifit(show_components=show_components, verbose=verbose,
                              debug=debug, use_lmfit=use_lmfit,
                              fittype=fittype,
                              guesses=guesses, annotate=annotate, **kwargs)
            else:
                raise ValueError("Guess and parinfo were somehow invalid.")
        else:
            raise ValueError("Can't fit with given fittype {0}:"
                             " it is not Registered as a fitter.".format(self.fittype))
        if save:
            self.savefit()

    def EQW(self, plot=False, plotcolor='g', fitted=True, continuum=None,
            components=False, annotate=False, alpha=0.5, loc='lower left',
            xmin=None, xmax=None, xunits=None, continuum_as_baseline=False,
            xunit='pixel', midpt_location='plot-center'):
        """
        Returns the equivalent width (integral of "baseline" or "continuum"
        minus the spectrum) over the selected range
        (the selected range defaults to self.xmin:self.xmax, so it may include
        multiple lines!)

        Parameters
        ----------
        plot : bool
            Plots a box indicating the EQW if plot==True (i.e., it will have a
            width equal to the equivalent width, and a height equal to the
            measured continuum)
        fitted : bool
            Use the fitted model?  If false, uses the data
        continuum : None or float
            Can specify a fixed continuum with this keyword, otherwise will use
            the fitted baseline.  WARNING: continuum=0 will still "work", but
            will give numerically invalid results.  Similarly, a negative continuum
            will work, but will yield results with questionable physical meaning.
        continuum_as_baseline : bool
            Replace the baseline with the specified continuum when computing
            the absorption depth of the line
        components : bool
            If your fit is multi-component, will attempt to acquire centroids
            for each component and print out individual EQWs
        xmin : float
        xmax : float
            The range over which to compute the EQW
        xunit : str
            The units of xmin/xmax
        midpt_location : 'fitted', 'plot-center'
            If 'plot' is set, this determines where the EQW will be drawn.  It
            can be the fitted centroid or the plot-center, i.e. (xmin+xmax)/2

        Returns
        -------
        Equivalent Width, or widths if components=True

        """
        if continuum is not None:
            # if continuum is specified, don't bother with checks
            if np.median(self.Spectrum.baseline.basespec) == 0:
                raise ValueError("Baseline / continuum is zero: equivalent width is undefined.")
            elif np.median(self.Spectrum.baseline.basespec) < 0:
                if mycfg.WARN: warn( "WARNING: Baseline / continuum is negative: equivalent width is poorly defined." )

        if xunits is not None and xunit=='pixel':
            # todo: deprecation warning
            xunit = xunits

        # determine range to use
        if xmin is None:
            xmin = self.xmin #self.Spectrum.xarr.x_to_pix(self.xmin)
        else:
            xmin = self.Spectrum.xarr.x_to_pix(xmin, xval_units=xunit)
        if xmax is None:
            xmax = self.xmax #self.Spectrum.xarr.x_to_pix(self.xmax)
        else:
            xmax = self.Spectrum.xarr.x_to_pix(xmax, xval_units=xunit)

        dx = np.abs(self.Spectrum.xarr[xmin:xmax].cdelt(approx=True).value)

        log.debug("xmin={0} xmax={1} dx={2} continuum={3}"
                  .format(xmin, xmax, dx, continuum))

        if components:
            centroids = self.fitter.analytic_centroids()
            integrals = self.fitter.component_integrals(self.Spectrum.xarr[xmin:xmax],dx=dx)
            eqw = []
            for cen,integ in zip(centroids,integrals):
                center_pix = self.Spectrum.xarr.x_to_pix(cen)
                if continuum is None:
                    continuum = self.Spectrum.baseline.basespec[center_pix]
                elif continuum_as_baseline:
                    integrals[-1] += -(self.Spectrum.baseline.basespec[xmin:xmax] - continuum).sum() * dx
                eqw.append(-integ / continuum)
            if plot:
                plot = False
                if mycfg.WARN:
                    warn("Cannot plot multiple Equivalent Widths")
        elif fitted:
            model = self.get_model(self.Spectrum.xarr[xmin:xmax],
                                   add_baseline=False)

            # EQW is positive for absorption lines
            # fitted components are assume to be continuum-subtracted
            integral = (-model).sum() * dx

            if continuum is None:
                # centroid in data units
                # (may fail if model has pos + neg values)
                center = (model*self.Spectrum.xarr[xmin:xmax]).sum()/model.sum()
                center_pix = self.Spectrum.xarr.x_to_pix(center)
                continuum = self.Spectrum.baseline.basespec[center_pix]
            elif continuum_as_baseline:
                integral += -(self.Spectrum.baseline.basespec[xmin:xmax] - continuum).sum() * dx

            eqw = integral / continuum
        else:
            if continuum_as_baseline:
                diffspec = (continuum - self.Spectrum.data)
            elif self.Spectrum.baseline.subtracted is False:
                diffspec = (self.Spectrum.baseline.basespec - self.Spectrum.data)
            else:
                diffspec = -self.Spectrum.data
            sumofspec = diffspec[xmin:xmax].sum() * dx
            if continuum is None:
                continuum = np.median(self.Spectrum.baseline.basespec)
            eqw = sumofspec / continuum
        if plot and self.Spectrum.plotter.axis:
            if midpt_location == 'plot-center':
                midpt_pixel = int(np.round((xmin+xmax)/2.0))
                midpt       = self.Spectrum.xarr[midpt_pixel].value
            elif midpt_location == 'fitted':
                try:
                    shifts = np.array([self.Spectrum.specfit.parinfo[x].value
                                       for x in self.Spectrum.specfit.parinfo.keys()
                                       if 'SHIFT' in x])
                except AttributeError:
                    raise AttributeError("Can only specify midpt_location="
                                         "fitted if there is a SHIFT parameter"
                                         "for the fitted model")
                # We choose to display the eqw fit at the center of the fitted
                # line set, closest to the passed window.
                # Note that this has the potential to show a eqw "rectangle"
                # centered on a fitted line other than the one measured for the
                # eqw call, if there are more than one fitted lines within the
                # window.
                midpt_pixel = int((xmin+xmax)/2)
                midval = self.Spectrum.xarr[midpt_pixel].value
                midpt_index = np.argmin(np.abs(shifts-midval))
                midpt = shifts[midpt_index]
                midpt_pixel = self.Spectrum.xarr.x_to_pix(midpt)
            else:
                raise ValueError("midpt_location must be 'plot-center' or "
                                 "fitted")
            if continuum_as_baseline:
                midpt_level = continuum
            else:
                midpt_level = self.Spectrum.baseline.basespec[midpt_pixel]
            log.debug("EQW plotting: midpt={0}, midpt_pixel={1}, "
                      "midpt_level={2}, eqw={3}".format(midpt, midpt_pixel,
                                                        midpt_level, eqw))
            self.EQW_plots.append(self.Spectrum.plotter.axis.fill_between(
                [midpt-eqw/2.0,midpt+eqw/2.0], [0,0],
                [midpt_level,midpt_level], color=plotcolor, alpha=alpha,
                label='EQW: %0.3g' % eqw))
            if annotate:
                self.Spectrum.plotter.axis.legend(
                        [(matplotlib.collections.CircleCollection([0],facecolors=[plotcolor],edgecolors=[plotcolor]))],
                        [('EQW: %0.3g' % eqw)],
                        markerscale=0.01, borderpad=0.1, handlelength=0.1,
                        handletextpad=0.1, loc=loc)
            if self.Spectrum.plotter.autorefresh:
                self.Spectrum.plotter.refresh()
        if hasattr(self.Spectrum,'header'):
            history.write_history(self.Spectrum.header, "EQW for %s: %s" %
                                  (self.fittype,eqw))
        return eqw

    def register_fitter(self,*args,**kwargs):
        """
        Register a model fitter
        """
        self.Registry.add_fitter(*args,**kwargs)

    register_fitter.__doc__ += Registry.add_fitter.__doc__

    def seterrspec(self, usestd=None, useresiduals=True):
        """
        Simple wrapper function to set the error spectrum; will either use the
        input spectrum or determine the error using the RMS of the residuals,
        depending on whether the residuals exist.
        """
        if (self.Spectrum.error is not None) and not usestd:
            if (self.Spectrum.error == 0).all():
                if self.residuals is not None and useresiduals:
                    residuals_std = self.residuals.std()
                    if residuals_std == 0:
                        self.errspec = np.ma.ones(self.spectofit.shape[0])
                        warnings.warn("Residuals have 0 standard deviation.  "
                                      "That's probably too good to be true.")
                    else:
                        self.errspec = np.ma.ones(self.spectofit.shape[0]) * residuals_std
                elif type(self.Spectrum.error) is np.ma.masked_array:
                    # force errspec to be a masked array of ones
                    self.errspec = self.Spectrum.error + 1
                else:
                    self.errspec = np.ma.copy(self.Spectrum.error) + 1
            else:
                # this is the default behavior if spectrum.error is set
                self.errspec = np.ma.copy(self.Spectrum.error)
        elif self.residuals is not None and useresiduals:
            self.errspec = np.ma.ones(self.spectofit.shape[0]) * self.residuals.std()
        else:
            self.errspec = np.ma.ones(self.spectofit.shape[0]) * self.spectofit.std()

        assert hasattr(self.errspec, 'mask')

    def setfitspec(self):
        """
        Set the spectrum that will be fit.  This is primarily to remove NANs
        from consideration: if you simply remove the data from both the X-axis
        and the Y-axis, it will not be considered for the fit, and a linear
        X-axis is not needed for fitting.

        However, it may be possible to do this using masked arrays instead of
        setting errors to be 1e10....
        """
        if self.Spectrum.data.sum() is np.ma.masked:
            self.spectofit = np.zeros_like(self.Spectrum.data)
            self.errspec = np.zeros_like(self.Spectrum.data)
            self._valid = False
            return
        # see https://github.com/numpy/numpy/issues/3474
        self.spectofit = np.ma.copy(self.Spectrum.data)
        if hasattr(self.Spectrum.data, 'mask') and hasattr(self.spectofit,
                                                           'mask'):
            assert np.all(self.Spectrum.data.mask == self.spectofit.mask)
        self._valid = True
        if hasattr(self.Spectrum,'baseline'):
            if ((not self.Spectrum.baseline.subtracted and
                 self.Spectrum.baseline.basespec is not None and
                 len(self.spectofit) == len(self.Spectrum.baseline.basespec))):
                self.spectofit -= self.Spectrum.baseline.basespec
        OKmask = np.isfinite(self.spectofit)

        with warnings.catch_warnings():
            # catch a specific np1.7 futurewarning relating to masks
            warnings.simplefilter("ignore")
            self.spectofit[~OKmask] = 0

        self.spectofit.mask = ~OKmask

        self.seterrspec()

        # the "OK" mask is just checking that the values are finite
        self.errspec[~OKmask] = 1e10
        self.errspec.mask = ~OKmask

        # if an includemask is set *and* there are some included values, "mask out" the rest
        # otherwise, if *all* data are excluded, we should assume that means the includemask
        # simply hasn't been initialized
        if self.includemask is not None and (self.includemask.shape == self.errspec.shape) and any(self.includemask):
            self.errspec[~self.includemask] = 1e10*self.errspec.max()

    @property
    def mask(self):
        """ Mask: True means "exclude" """
        if ((hasattr(self.spectofit, 'mask') and
             self.spectofit.shape==self.spectofit.mask.shape)):
            mask = self.spectofit.mask
        else:
            mask = np.zeros_like(self.spectofit, dtype='bool')

        return mask

    @property
    def dof(self):
        """ degrees of freedom in fit """
        if not hasattr(self, 'npix_fitted'):
            raise AttributeError('No fit has been run, so npix_fitted is not '
                                 'defined and dof cannot be computed.')
        return (self.npix_fitted - self.vheight - self.npeaks *
                self.Registry.npars[self.fittype] + np.sum(self.parinfo.fixed) +
                np.sum([x != '' for x in self.parinfo.tied]))

        #self.dof  = self.includemask.sum()-self.npeaks*self.Registry.npars[self.fittype]-vheight+np.sum(self.parinfo.fixed)


    @property
    def mask_sliced(self):
        """ Sliced (subset) Mask: True means "exclude" """
        return self.mask[self.xmin:self.xmax]

    def multifit(self, fittype=None, renormalize='auto', annotate=None,
                 show_components=None, verbose=True, color=None,
                 guesses=None, parinfo=None, reset_fitspec=True,
                 use_window_limits=None, use_lmfit=False, plot=True, **kwargs):
        """
        Fit multiple gaussians (or other profiles)

        Parameters
        ----------
        fittype : str
            What function will be fit?  fittype must have been Registryed in the
            peakbgfitters dict.  Uses default ('gaussian') if not specified
        renormalize : 'auto' or bool
            if 'auto' or True, will attempt to rescale small data (<1e-9) to be
            closer to 1 (scales by the median) so that the fit converges better
        parinfo : `~parinfo` structure
            Guess structure; supercedes ``guesses``
        guesses : list or 'moments'
            Either a list of guesses matching the number of parameters * the
            number of peaks for the model, or 'moments' to fit a single
            spectrum with the moments as guesses

        """
        if reset_fitspec:
            self.setfitspec()
        if not self._valid:
            raise ValueError("Data are invalid; cannot be fit.")
        #if self.fitkwargs.has_key('negamp'): self.fitkwargs.pop('negamp') # We now do this in gaussfitter.py
        if fittype is not None:
            self.fittype = fittype
        bad_kws = ['fittype','plot']
        for kw in bad_kws:
            if kw in self.fitkwargs:
                del self.fitkwargs[kw]

        if guesses is not None and parinfo is not None:
            raise ValueError("Both guesses and parinfo were specified, "
                             "but only one of these is allowed.")

        if guesses is None:
            if parinfo is not None:
                guesses = list(parinfo.values)
            else:
                guesses = list(self.guesses)
        elif isinstance(guesses, string_types) and guesses in ('moment', 'moments'):
            guesses = self.moments(vheight=False, **kwargs)
        else:
            guesses = list(guesses) # needs to be mutable, but needs to be a copy!!


        if len(guesses) < self.Registry.npars[self.fittype]:
            raise ValueError("Too few parameters input.  Need at least %i for %s models" % (self.Registry.npars[self.fittype],self.fittype))

        self.npeaks = len(guesses)/self.Registry.npars[self.fittype]
        self.fitter = self.Registry.multifitters[self.fittype]
        self.vheight = False
        if self.fitter.vheight:
            # Need to reset the parinfo if vheight has previously been set,
            # otherwise npars will disagree, which causes problems if
            # renormalization happens
            self.fitter.vheight = False
            self.fitter.npeaks = self.npeaks
            self.fitter._make_parinfo(npeaks=self.npeaks)

        # add kwargs to fitkwargs
        self.fitkwargs.update(kwargs)
        if 'renormalize' in self.fitkwargs:
            del self.fitkwargs['renormalize']

        # if parinfo was specified, we use it and ignore guesses otherwise, we
        # make a parinfo so we can test 'scaleable' below
        if parinfo is not None:
            pinf_for_scaling = parinfo
        else:
            pinf_for_scaling, _ = self.fitter._make_parinfo(parvalues=guesses,
                                                            npeaks=self.npeaks,
                                                            **self.fitkwargs)

        scalefactor = 1.0
        if renormalize in ('auto',True):
            datarange = np.nanmax(self.spectofit[self.xmin:self.xmax]) - np.nanmin(self.spectofit[self.xmin:self.xmax])
            if abs(datarange) < 1e-9:
                scalefactor = np.nanmedian(np.abs(self.spectofit))
                if not np.isfinite(scalefactor):
                    raise ValueError("non-finite scalefactor = {0} encountered.".format(scalefactor))
                elif scalefactor == 0:
                    raise ValueError("scalefactor = {0} encountered, which will result "
                                     "in divide-by-zero errors".format(scalefactor))
                log.info("Renormalizing data by factor %e to improve fitting procedure"
                         % scalefactor)
                self.spectofit /= scalefactor
                self.errspec /= scalefactor

                # this error should be unreachable, but is included as a sanity check
                if self.fitter.npeaks * self.fitter.npars != len(pinf_for_scaling):
                    raise ValueError("Length of parinfo doesn't agree with "
                                     " npeaks * npars = {0}"
                                     .format(self.fitter.npeaks *
                                             self.fitter.npars))
                if len(guesses) != len(pinf_for_scaling):
                    raise ValueError("Length of parinfo doens't match length of guesses")

                # zip guesses with parinfo: truncates parinfo if len(parinfo) > len(guesses)
                # actually not sure how/when/if this should happen; this might be a bad hack
                # revisit with tests!!
                for jj,(guess,par) in enumerate(zip(guesses,pinf_for_scaling)):
                    if par.scaleable:
                        guesses[jj] /= scalefactor
                        # if parinfo was passed in, this will change it
                        # if it was not, it will change only the placeholder
                        # (becuase we are passing by reference above)
                        par.value /= scalefactor
                        par.limits = [lim / scalefactor for lim in par.limits]

                log.debug("Rescaled guesses to {0}".format(guesses))

        # all fit data must be float64, otherwise the optimizers may take steps
        # less than the precision of the data and get stuck
        xtofit = self.Spectrum.xarr[self.xmin:self.xmax][~self.mask_sliced].astype('float64')
        spectofit = self.spectofit[self.xmin:self.xmax][~self.mask_sliced].astype('float64')
        err = self.errspec[self.xmin:self.xmax][~self.mask_sliced].astype('float64')

        if hasattr(xtofit, 'value') and not hasattr(xtofit, 'x_to_coord'):
            xtofit = SpectroscopicAxis(xtofit)

        if np.all(err == 0):
            raise ValueError("Errors are all zero.  This should not occur and "
                             "is a bug. (if you set the errors to all zero, "
                             "they should be overridden and set to 1)")

        if parinfo is not None:
            self._validate_parinfo(parinfo, mode='fix')
        else:
            pinf, _ = self.fitter._make_parinfo(parvalues=guesses,
                                                npeaks=self.npeaks,
                                                **self.fitkwargs)
            new_guesses = self._validate_parinfo(pinf, 'guesses')
            if any((x!=y) for x,y in zip(guesses, new_guesses)):
                warn("Guesses have been changed from {0} to {1}"
                     .format(guesses, new_guesses))
            guesses = new_guesses

        mpp,model,mpperr,chi2 = self.fitter(xtofit, spectofit, err=err,
                                            npeaks=self.npeaks,
                                            parinfo=parinfo, # the user MUST be allowed to override parinfo.
                                            params=guesses,
                                            use_lmfit=use_lmfit,
                                            **self.fitkwargs)

        any_out_of_range = self._validate_parinfo(self.fitter.parinfo, mode='check')

        if any(any_out_of_range):
            warn("The fitter returned values that are outside the "
                 "parameter limits.  DEBUG INFO: {0}".format(any_out_of_range))

        self.spectofit *= scalefactor
        self.errspec *= scalefactor

        if hasattr(self.fitter.mp,'status'):
            self.mpfit_status = models.mpfit_messages[self.fitter.mp.status]

        if model is None:
            raise ValueError("Model was not set by fitter.  Examine your fitter.")
        self.chi2 = chi2
        self.model = model * scalefactor
        self.parinfo = self.fitter.parinfo

        # rescale any scaleable parameters
        for par in self.parinfo:
            if par.scaleable:
                par.value *= scalefactor
                if par.error is not None:
                    par.error *= scalefactor
                if par.limits is not None:
                    par.limits = [lim*scalefactor for lim in par.limits]

        self.modelpars = self.parinfo.values
        self.modelerrs = self.parinfo.errors
        self.residuals = spectofit - self.model
        if self.Spectrum.plotter.axis is not None and plot:
            if color is not None:
                kwargs.update({'composite_fit_color':color})
            self.plot_fit(annotate=annotate,
                          show_components=show_components,
                          use_window_limits=use_window_limits,
                          **kwargs)

        # Re-organize modelerrs so that any parameters that are tied to others inherit the errors of the params they are tied to
        if 'tied' in self.fitkwargs:
            for ii, element in enumerate(self.fitkwargs['tied']):
                if not element.strip():
                    continue

                if '[' in element and ']' in element:
                    i1 = element.index('[') + 1
                    i2 = element.index(']')
                    loc = int(element[i1:i2])
                else: # assume lmfit version
                    varnames = re.compile('([a-zA-Z][a-zA-Z_0-9]*)').search(element).groups()
                    if not varnames:
                        continue
                    elif len(varnames) > 1:
                        warn("The 'tied' parameter {0} is not simple enough for error propagation".format(element))
                        continue
                    else:
                        varname = varnames[0]
                        loc = self.parinfo.names.index(varname)

                self.modelerrs[ii] = self.modelerrs[loc]

        # make sure the full model is populated
        self._full_model()

        # calculate the number of pixels included in the fit.  This should
        # *only* be done when fitting, not when selecting data.
        # (see self.dof)
        self.npix_fitted = self.includemask.sum() - self.mask.sum()

        self.history_fitpars()


    def refit(self, use_lmfit=False):
        """ Redo a fit using the current parinfo as input """
        return self.multifit(parinfo=self.parinfo, use_lmfit=use_lmfit,
                             reset_fitspec=False)

    def history_fitpars(self):
        if hasattr(self.Spectrum,'header'):
            history.write_history(self.Spectrum.header, "SPECFIT: Fitted "
                                  "profile of type %s" % (self.fittype))
            history.write_history(self.Spectrum.header, "Chi^2: %g  DOF: %i" %
                                  (self.chi2, self.dof))
            for par in self.parinfo:
                history.write_history(self.Spectrum.header, str(par))

    def peakbgfit(self, usemoments=True, annotate=None, vheight=True, height=0,
                  negamp=None, fittype=None, renormalize='auto', color=None,
                  use_lmfit=False, show_components=None, debug=False,
                  use_window_limits=True, guesses=None,
                  nsigcut_moments=None, plot=True, parinfo=None, **kwargs):
        """
        Fit a single peak (plus a background)

        Parameters
        ----------
        usemoments : bool
            The initial guess will be set by the fitter's 'moments' function
            (this overrides 'guesses')
        annotate : bool
            Make a legend?
        vheight : bool
            Fit a (constant) background as well as a peak?
        height : float
            initial guess for background
        negamp : bool
            If True, assumes amplitude is negative.  If False, assumes positive.  If
            None, can be either.
        fittype : bool
            What function will be fit?  fittype must have been Registryed in the
            peakbgfitters dict
        renormalize : 'auto' or bool
            if 'auto' or True, will attempt to rescale small data (<1e-9) to be
            closer to 1 (scales by the median) so that the fit converges better
        nsigcut_moments : bool
            pass to moment guesser; can do a sigma cut for moment guessing


        """
        self.npeaks = 1
        self.auto = True
        self.setfitspec()

        if fittype is not None:
            self.fittype=fittype
        NP = self.Registry.peakbgfitters[self.fittype].default_npars

        if guesses is not None:
            log.debug("Using user-specified guesses.")
            self.guesses = guesses
            if len(guesses) != NP + vheight:
                raise ValueError("Invalid guesses specified for single-fitter."
                                 "Expected {0}, got {1}.  Perhaps you should "
                                 "use the multifitter (multifit=True)?"
                                 .format(NP+vheight, len(guesses)))

        elif usemoments: # this can be done within gaussfit but I want to save them
            # use this INDEPENDENT of fittype for now (voigt and gauss get same guesses)
            log.debug("Using moment-based guesses.")
            moments_f = self.Registry.peakbgfitters[self.fittype].moments
            self.guesses = moments_f(self.Spectrum.xarr[self.xmin:self.xmax],
                                     self.spectofit[self.xmin:self.xmax],
                                     vheight=vheight,
                                     negamp=negamp,
                                     nsigcut=nsigcut_moments,
                                     **kwargs)
        else:
            if negamp:
                self.guesses = [height,-1,0,1]
            else:
                self.guesses = [height,1,0,1]

        # If we're fitting anything but a simple Gaussian, we need the length
        # of guesses to be right so we pad with appended zeros
        # BUT, if the guesses from the moments have the right number of
        # parameters, we don't need to do this.
        if NP > len(self.guesses):
            for ii in xrange(len(self.guesses),NP):
                self.guesses += [0.0]

        self.fitter = self.Registry.peakbgfitters[self.fittype]

        log.debug("n(guesses): %s  Guesses: %s  vheight: %s " %
                  (len(self.guesses),self.guesses,vheight))

        scalefactor = 1.0
        if renormalize in ('auto',True):
            datarange = self.spectofit[self.xmin:self.xmax].max() - self.spectofit[self.xmin:self.xmax].min()
            if abs(datarange) < 1e-9:
                scalefactor = np.median(np.abs(self.spectofit))
                log.info("Renormalizing data by factor %e to improve fitting procedure"
                         % scalefactor)
                self.spectofit /= scalefactor
                self.errspec /= scalefactor
                self.guesses[0] /= scalefactor
                if vheight:
                    self.guesses[1] /= scalefactor

        log.debug("Guesses before fit: {0}".format(self.guesses))

        if 'debug' in self.fitkwargs:
            debug = self.fitkwargs['debug']
            del self.fitkwargs['debug']

        mpp,model,mpperr,chi2 = self.fitter(
            self.Spectrum.xarr[self.xmin:self.xmax],
            self.spectofit[self.xmin:self.xmax],
            err=self.errspec[self.xmin:self.xmax], vheight=vheight,
            params=self.guesses, parinfo=parinfo, debug=debug,
            use_lmfit=use_lmfit, **self.fitkwargs)
        log.debug("1. Guesses, fits after: {0}, {1}".format(self.guesses, mpp))

        self.spectofit *= scalefactor
        self.errspec *= scalefactor

        if hasattr(self.fitter.mp,'status'):
            self.mpfit_status = models.mpfit_messages[self.fitter.mp.status]
        self.parinfo = self.fitter.parinfo

        if model is None:
            raise ValueError("Model was not set by fitter.  Examine your fitter.")
        self.chi2 = chi2
        self.vheight=vheight
        if vheight:
            self.Spectrum.baseline.order = 0
            self.Spectrum.baseline.baselinepars = [mpp[0]*scalefactor] # first item in list form
            self.Spectrum.baseline.basespec = self.Spectrum.data*0 + mpp[0]*scalefactor
            self.model = model*scalefactor - mpp[0]*scalefactor
            # I removed this recently for some reason, but more code depends on it being in place
            # Need to figure out *WHY* anything would want an extra parameter
            if len(mpp) == self.fitter.npars+1:
                mpp = mpp[1:]
        else:
            self.model = model*scalefactor
        self.residuals = self.spectofit[self.xmin:self.xmax] - self.model*scalefactor
        self.modelpars = mpp
        self.modelerrs = mpperr

        # rescale any scaleable parameters
        for par in self.parinfo:
            if par.scaleable:
                par.value = par.value * scalefactor
                if par.error is not None:
                    par.error = par.error * scalefactor

        if self.Spectrum.plotter.axis is not None and plot:
            if color is not None:
                kwargs.update({'composite_fit_color':color})
            self.plot_fit(annotate=annotate,
                          use_window_limits=use_window_limits,
                          show_components=show_components,
                          **kwargs)

        # make sure the full model is populated
        self._full_model(debug=debug)

        self.npix_fitted = self.includemask.sum() - self.mask.sum()

        log.debug("2. Guesses, fits after vheight removal: {0},{1}"
                  .format(self.guesses, mpp))
        self.history_fitpars()

    def _full_model(self, debug=False, **kwargs):
        """
        Compute the model for the whole spectrum
        """
        self.fullmodel = self.get_full_model(debug=debug,**kwargs)
        self.fullresiduals = self.Spectrum.data - self.fullmodel

    def get_full_model(self, debug=False,**kwargs):
        """ compute the model over the full axis """
        return self.get_model(self.Spectrum.xarr, debug=debug,**kwargs)

    def get_model(self, xarr, pars=None, debug=False, add_baseline=None):
        """ Compute the model over a given axis """
        if pars is None:
            return self.get_model_frompars(xarr=xarr, pars=self.parinfo,
                    add_baseline=add_baseline, debug=debug)
        else:
            return self.get_model_frompars(xarr=xarr, pars=pars,
                    add_baseline=add_baseline, debug=debug)

    def get_model_frompars(self, xarr, pars, debug=False, add_baseline=None):
        """ Compute the model over a given axis """
        if ((add_baseline is None and (self.Spectrum.baseline.subtracted or self.vheight))
                or add_baseline is False):
            return self.fitter.n_modelfunc(pars,**self.fitter.modelfunc_kwargs)(xarr)
        else:
            return (self.fitter.n_modelfunc(pars,
                                            **self.fitter.modelfunc_kwargs)(xarr)
                    + self.Spectrum.baseline.get_model(np.arange(xarr.size)))

    def plot_model(self, pars, offset=0.0, annotate=False, clear=False, **kwargs):
        """
        Plot a model from specified input parameters
        (see plot_fit for kwarg specification)

        annotate is set to "false" because arbitrary annotations are not yet implemented
        """

        # really, plot_fit should be thin on top of plot_model, but that's
        # not how I wrote it, so it will have to wait for a refactor

        if clear: self.clear()

        return self.plot_fit(pars=pars, offset=offset, annotate=False, **kwargs)


    #def assess_npeaks(self):
    #    """
    #    Attempt to determine whether any of the peaks are unnecessary
    #    """
    #    if self.npeaks <= 1:
    #        return
    #    npars = self.fitter.npars
    #    perpeakpars = [self.parinfo.values[ii*npars:(ii+1)*npars] for ii in
    #                   range(self.npeaks)]
    #    parsets = [((x[0][0],x[1][0]),x[0][1]+x[1][1]) for x in
    #               itertools.combinations(perpeakpars, self.npeaks-1)]
    #    parsets = [x
    #               for y in itertools.combinations(perpeakpars, self.npeaks-1)
    #               for x in y]

    #    chi2_without = [(self.spectofit[self.xmin:self.xmax] -
    #                     self.get_model_frompars(self.xarr, self.pars[ii*npars:

    def plot_fit(self, xarr=None, annotate=None, show_components=None,
                 composite_fit_color='red',  lw=0.5,
                 composite_lw=0.75, pars=None, offset=None,
                 use_window_limits=None, show_hyperfine_components=None,
                 plotkwargs={}, **kwargs):
        """
        Plot the fit.  Must have fitted something before calling this!

        It will be automatically called whenever a spectrum is fit (assuming an
        axis for plotting exists)

        kwargs are passed to the fitter's components attribute

        Parameters
        ----------
        xarr : None
            If none, will use the spectrum's xarr.  Otherwise, plot the
            specified xarr.  This is useful if you want to plot a well-sampled
            model when the input spectrum is undersampled
        annotate : None or bool
            Annotate the plot?  If not specified, defaults to self.autoannotate
        show_components : None or bool
        show_hyperfine_components : None or bool
            Show the individual gaussian components overlaid on the composite fit
        use_window_limits : None or bool
            If False, will reset the window to include the whole spectrum.  If
            True, leaves the window as is.  Defaults to self.use_window_limits
            if None.
        pars : parinfo
            A parinfo structure or list of model parameters.  If none, uses
            best-fit
        offset : None or float
            Y-offset. If none, uses the default self.Spectrum.plotter offset, otherwise,
            uses the specified float.
        """
        #if self.Spectrum.baseline.subtracted is False and self.Spectrum.baseline.basespec is not None:
        #    # don't display baseline if it's included in the fit
        #    plot_offset = self.Spectrum.plotter.offset+(self.Spectrum.baseline.basespec * (~self.vheight))
        #else:
        if offset is None:
            plot_offset = self.Spectrum.plotter.offset
        else:
            plot_offset = offset

        if xarr is None:
            xarr = self.Spectrum.xarr

        if pars is not None:
            model = self.get_model_frompars(xarr, pars)
        else:
            self._full_model()
            model = self.fullmodel

        self.modelplot += self.Spectrum.plotter.axis.plot(xarr,
                                                          model + plot_offset,
                                                          color=composite_fit_color,
                                                          linewidth=lw,
                                                          **plotkwargs)

        # Plot components
        if show_components or show_hyperfine_components:
            self.plot_components(xarr=xarr,
                                 show_hyperfine_components=show_hyperfine_components,
                                 pars=pars, plotkwargs=plotkwargs)

        uwl = use_window_limits if use_window_limits is not None else self.use_window_limits
        # plotter kwargs are kwargs for the Spectrum.plotter,
        # whereas plotkwargs are for the matplotlib plot command
        plotterkwargs = {}
        plotterkwargs.update(self.Spectrum.plotter.plotkwargs)
        plotterkwargs['use_window_limits'] = uwl
        self.Spectrum.plotter.reset_limits(**plotterkwargs)
        if self.Spectrum.plotter.autorefresh:
            self.Spectrum.plotter.refresh()

        if annotate or ((annotate is None) and self.autoannotate):
            self.annotate()
            if self.vheight: self.Spectrum.baseline.annotate()

    def plot_components(self, xarr=None, show_hyperfine_components=None,
            component_yoffset=0.0, component_lw=0.75, pars=None,
            component_fit_color='blue', component_kwargs={},
            add_baseline=False, plotkwargs={}, **kwargs):
        """
        Overplot the individual components of a fit

        Parameters
        ----------
        xarr : None
            If none, will use the spectrum's xarr.  Otherwise, plot the
            specified xarr.  This is useful if you want to plot a well-sampled
            model when the input spectrum is undersampled
        show_hyperfine_components : None | bool
            Keyword argument to pass to component codes; determines whether to return
            individual (e.g., hyperfine) components of a composite model
        component_yoffset : float
            Vertical (y-direction) offset to add to the components when plotting
        component_lw : float
            Line width of component lines
        component_fitcolor : color
            Color of component lines
        component_kwargs : dict
            Keyword arguments to pass to the fitter.components method
        add_baseline : bool
            Add the fit to the components before plotting.  Makes sense to use
            if self.Spectrum.baseline.subtracted == False
        pars : parinfo
            A parinfo structure or list of model parameters.  If none, uses
            best-fit
        """

        plot_offset = self.Spectrum.plotter.offset

        if xarr is None:
            xarr = self.Spectrum.xarr

        if show_hyperfine_components is not None:
            component_kwargs['return_hyperfine_components'] = show_hyperfine_components
        self._component_kwargs = component_kwargs

        if pars is None:
            pars = self.modelpars

        self.modelcomponents = self.fitter.components(xarr=xarr, pars=pars, **component_kwargs)

        yoffset = plot_offset + component_yoffset
        if add_baseline:
            yoffset += self.Spectrum.baseline.basespec

        for data in self.modelcomponents:
            # can have multidimensional components
            if len(data.shape) > 1:
                for d in data:
                    self._plotted_components += self.Spectrum.plotter.axis.plot(xarr,
                        d + yoffset,
                        color=component_fit_color, linewidth=component_lw, **plotkwargs)
            else:
                self._plotted_components += self.Spectrum.plotter.axis.plot(xarr,
                    data + yoffset,
                    color=component_fit_color, linewidth=component_lw, **plotkwargs)


    def fullsizemodel(self):
        """
        If the model was fit to a sub-region of the spectrum, expand it (with
        zeros wherever the model was not defined) to fill the spectrum.

        Examples
        --------
        >>> noise = np.random.randn(100)
        >>> xarr = np.linspace(-50,50,100)
        >>> signal = np.exp(-(xarr-5)**2/(2*3.**2))
        >>> sp = pyspeckit.Spectrum(data=noise + signal, xarr=xarr, xarrkwargs={'units':'km/s'})
        >>> sp.specfit(xmin=-25,xmax=25)
        >>> sp.specfit.model.shape
        (48,)
        >>> sp.specfit.fullsizemodel()
        >>> sp.specfit.model.shape
        (100,)
        """

        if self.model.shape != self.Spectrum.data.shape:
            temp = np.zeros(self.Spectrum.data.shape)
            temp[self.xmin:self.xmax] = self.model
            self.model = temp
            self.residuals = self.spectofit - self.model
            self.selectregion(reset=True)

    def plotresiduals(self, fig=2, axis=None, clear=True, color='k',
                      linewidth=0.5, drawstyle='steps-mid', yoffset=0.0,
                      label=True, pars=None, zeroline=None,
                      set_limits=True, **kwargs):
        """
        Plot residuals of the fit.  Specify a figure or
        axis; defaults to figure(2).

        Parameters
        ----------
        fig : int
            Figure number.  Overridden by axis
        axis : axis
            The axis to plot on
        pars : None or parlist
            If set, the residuals will be computed for the input parameters
        zeroline : bool or None
            Plot the "zero" line through the center of the residuals.  If None,
            defaults to "True if yoffset!=0, False otherwise"


        kwargs are passed to matplotlib plot
        """
        self._full_model(pars=pars)
        if axis is None:
            if isinstance(fig,int):
                fig=matplotlib.pyplot.figure(fig)
            self.residualaxis = matplotlib.pyplot.gca()
            if clear:
                self.residualaxis.clear()
        else:
            self.residualaxis = axis
            if clear:
                self.residualaxis.clear()
        self.residualplot = self.residualaxis.plot(self.Spectrum.xarr,
                                                   self.fullresiduals+yoffset,
                                                   drawstyle=drawstyle,
                                                   linewidth=linewidth,
                                                   color=color, **kwargs)
        if zeroline or (zeroline is None and yoffset != 0):
            self.residualplot += self.residualaxis.plot(self.Spectrum.xarr,
                                                        (np.zeros_like(self.Spectrum.xarr.value)+yoffset),
                                                        linestyle='--',
                                                        color='k',
                                                        alpha=0.5)
        if set_limits:
            if ((self.Spectrum.plotter.xmin is not None) and
                (self.Spectrum.plotter.xmax is not None)):
                self.residualaxis.set_xlim(self.Spectrum.plotter.xmin.value,
                                           self.Spectrum.plotter.xmax.value)
            if ((self.Spectrum.plotter.ymin is not None) and
                (self.Spectrum.plotter.ymax is not None)):
                self.residualaxis.set_ylim(self.Spectrum.plotter.ymin,
                                           self.Spectrum.plotter.ymax)
        if label:
            self.residualaxis.set_xlabel(self.Spectrum.plotter.xlabel)
            self.residualaxis.set_ylabel(self.Spectrum.plotter.ylabel)
            self.residualaxis.set_title("Residuals")
        if self.Spectrum.plotter.autorefresh:
            self.residualaxis.figure.canvas.draw()

    def annotate(self,loc='upper right',labelspacing=0.25, markerscale=0.01,
                 borderpad=0.1, handlelength=0.1, handletextpad=0.1,
                 fontsize=10,
                 frameon=False, chi2=None, optimal_chi2_kwargs={}, **kwargs):
        """
        Add a legend to the plot showing the fitted parameters

        _clearlegend() will remove the legend

        chi2 : {True or 'reduced' or 'optimal' or 'allthree'}

        kwargs passed to legend
        """
        self._clearlegend()
        pl = matplotlib.collections.CircleCollection([0],edgecolors=['k'])
        if hasattr(self.fitter,'annotations'):
            self._annotation_labels = self.fitter.annotations()
        else:
            raise Exception("Fitter %s has no annotations." % self.fitter)

        #xtypename = units.unit_type_dict[self.Spectrum.xarr.xtype]
        xcharconv = units.SmartCaseNoSpaceDict({u.Hz.physical_type:'\\nu',
                                                u.m.physical_type:'\\lambda',
                                                (u.km/u.s).physical_type:'v',
                                                'pixels':'x',
                                                u.dimensionless_unscaled:'x',
                                                'dimensionless':'x',
                                               })
        try:
            xchar = xcharconv[self.Spectrum.xarr.unit.physical_type]
        except AttributeError:
            unit_key = self.Spectrum.xarr.unit
            xchar = xcharconv[u.Unit(unit_key).physical_type]

        self._annotation_labels = [L.replace('x',xchar) if L[1]=='x' else L for
                                   L in self._annotation_labels]

        if chi2 is not None:
            chi2n_label = '$\\chi^2/\\nu = %0.2g$' % (self.chi2/self.dof)
            chi2opt_label = '$\\chi^2_o/\\nu = %0.2g$' % self.optimal_chi2(**optimal_chi2_kwargs)
            chi2_label = '$\\chi^2 = %0.2g$' % self.chi2
            if chi2 == 'allthree':
                self._annotation_labels.append("\n".join([chi2n_label,
                                                          chi2_label,
                                                          chi2opt_label]))
            elif chi2 == 'reduced':
                self._annotation_labels.append(chi2n_label)
            elif chi2 == 'optimal':
                self._annotation_labels.append(chi2opt_label)
            else:
                self._annotation_labels.append(chi2_label)

        if self.Spectrum.plotter.axis:
            try:
                self.fitleg = self.Spectrum.plotter.axis.legend(
                    tuple([pl]*len(self._annotation_labels)),
                    self._annotation_labels, loc=loc, markerscale=markerscale,
                    borderpad=borderpad, handlelength=handlelength,
                    handletextpad=handletextpad, labelspacing=labelspacing,
                    frameon=frameon, fontsize=fontsize, **kwargs)
                self.Spectrum.plotter.axis.add_artist(self.fitleg)
            except TypeError as ex:
                print("Error {0} was raised, which may indicate an outdated mpl version".format(ex))
            try:
                self.fitleg.set_draggable(True)
            except AttributeError:
                # wrong version and/or non-interactive backend
                pass
            if self.Spectrum.plotter.autorefresh:
                self.Spectrum.plotter.refresh()

    def print_fit(self, print_baseline=True, **kwargs):
        """
        Print the best-fit parameters to the command line
        """

        if self.Spectrum.baseline.baselinepars is not None and print_baseline:
            print("Baseline: " + " + ".join(["%12g x^%i" % (x,i) for i,x in enumerate(self.Spectrum.baseline.baselinepars[::-1])]))

        for i,p in enumerate(self.parinfo):
            print("%15s: %12g +/- %12g" % (p['parname'],p['value'],p['error']))

    def clear(self, legend=True, components=True):
        """
        Remove the fitted model from the plot

        Also removes the legend by default
        """
        if self.Spectrum.plotter.axis is not None:
            for p in self.modelplot:
                p.set_visible(False)
            if legend:
                self._clearlegend()
            if components:
                self._clearcomponents()
            if self.Spectrum.plotter.autorefresh:
                self.Spectrum.plotter.refresh()

        # Empty the modelplot array to free memory
        self.modelplot = []

        # remove residuals from self if they're there.
        if hasattr(self,'residualplot'):
            for L in self.residualplot:
                if L in self.Spectrum.plotter.axis.lines:
                    if hasattr(self.Spectrum.plotter.axis.lines, 'remove'):
                        self.Spectrum.plotter.axis.lines.remove(L)
                    else:
                        L.remove()

    def _clearcomponents(self):
        for pc in self._plotted_components:
            pc.set_visible(False)
            if pc in self.Spectrum.plotter.axis.lines:
                if hasattr(self.Spectrum.plotter.axis.lines, 'remove'):
                    self.Spectrum.plotter.axis.lines.remove(pc)
                else:
                    pc.remove()
        if self.Spectrum.plotter.autorefresh:
            self.Spectrum.plotter.refresh()

        # Empty the plotted components array to free memory
        self._plotted_components = []

    def _clearlegend(self):
        """
        Remove the legend from the plot window
        """
        axis = self.Spectrum.plotter.axis
        if axis and axis.legend_ == self.fitleg:
            axis.legend_ = None
        if axis and self.fitleg is not None:
            # don't remove fitleg unless it's in the current axis
            # self.fitleg.set_visible(False)
            if self.fitleg in axis.artists:
                try:
                    axis.artists.remove(self.fitleg)
                except AttributeError:
                    try:
                        self.fitleg.remove()
                    except Exception as ex:
                        pass
        if self.Spectrum.plotter.autorefresh:
            self.Spectrum.plotter.refresh()


    def savefit(self):
        """
        Save the fit parameters from a Gaussian fit to the FITS header

        .. todo::
            THESE SHOULD BE WRITTEN FOR EACH TYPE OF MODEL TO BE FIT
        """
        if self.modelpars is not None and hasattr(self.Spectrum,'header'):
            for ii,p in enumerate(self.modelpars):

                try:
                    if ii % 3 == 0:
                        self.Spectrum.header['AMP%1i' % (ii/3)] = (p,"Gaussian best fit amplitude #%i" % (ii/3))
                    elif ii % 3 == 1:
                        self.Spectrum.header['CEN%1i' % (ii/3)] = (p,"Gaussian best fit center #%i" % (ii/3))
                    elif ii % 3 == 2:
                        self.Spectrum.header['WID%1i' % (ii/3)] = (p,"Gaussian best fit width #%i" % (ii/3))
                except ValueError as ex:
                    log.info("Failed to save fit to header: {0}".format(ex))

    def downsample(self,factor):
        """
        Downsample the model spectrum (and the spectofit spectra)
        This should only be done when Spectrum.smooth is called
        """
        if self.model is not None:
            self.model = self.model[::factor]
        if self.residuals is not None:
            self.residuals = self.residuals[::factor]
        self.spectofit = self.spectofit[::factor]
        self.errspec = self.errspec[::factor]
        self.includemask = self.includemask[::factor]

    def crop(self,x1pix,x2pix):
        """
        When spectrum.crop is called, this must be too
        """
        if self.model is not None:
            self.model = self.model[x1pix:x2pix]
        if hasattr(self,'fullmodel'):
            self.fullmodel = self.fullmodel[x1pix:x2pix]
        self.includemask = self.includemask[x1pix:x2pix]
        self.setfitspec()

    def integral(self, analytic=False, direct=False, threshold='auto',
            integration_limits=None, integration_limit_units='pixels',
            return_error=False, **kwargs):
        """
        Return the integral of the fitted spectrum

        Parameters
        ----------
        analytic : bool
            Return the analytic integral of the fitted function?
            .. WARNING:: This approach is only implemented for some models
            .. todo:: Implement error propagation for this approach
        direct : bool
            Return the integral of the *spectrum* (as opposed to the *fit*)
            over a range defined by the `integration_limits` if specified or
            `threshold` otherwise
        threshold : 'auto' or 'error' or float
            Determines what data to be included in the integral based off of where
            the model is greater than this number
            If 'auto', the threshold will be set to peak_fraction * the peak
            model value.
            If 'error', uses the error spectrum as the threshold
            See `self.get_model_xlimits` for details
        integration_limits : None or 2-tuple
            Manually specify the limits in `integration_limit_units` units
        return_error : bool
            Return the error on the integral if set.
            The error computed by
            sigma = sqrt(sum(sigma_i^2)) * dx
        kwargs :
            passed to `self.fitter.integral` if ``not(direct)``

        Returns
        -------
        np.scalar or np.ndarray with the integral or integral & error

        """

        if analytic:
            return self.fitter.analytic_integral(modelpars=self.parinfo.values)

        xmin,xmax = self.get_model_xlimits(units='pixels', threshold=threshold)
        if integration_limits is None:
            integration_limits = [xmin,xmax]
        integration_limits = [
                    self.Spectrum.xarr.x_to_pix(x,xval_units=integration_limit_units)
                    for x in integration_limits]

        if xmax - xmin > 1: # can only get cdelt if there's more than 1 pixel
            dx = self.Spectrum.xarr[xmin:xmax].cdelt().value
        else:
            dx = None
        if dx is None:
            #dx = np.abs(np.concatenate([np.diff(self.Spectrum.xarr),[0]]))
            #warn("Irregular X-axis.  The last pixel is ignored.")
            self.Spectrum.xarr.make_dxarr()
            dx = self.Spectrum.xarr.dxarr.value
        else:
            # shouldn't shape be a 'property'?
            dx = np.repeat(np.abs(dx), self.Spectrum.shape)


        if direct:
            integrand = self.Spectrum.data[xmin:xmax]
            if not self.Spectrum.baseline.subtracted:
                integrand -= self.Spectrum.baseline.basespec[xmin:xmax]

            integ = (integrand * dx[xmin:xmax]).sum()
            if return_error:
                # compute error assuming a "known mean" (not a sample mean).  If sample mean, multiply
                # by sqrt(len(dx)/(len(dx)-1))  (which should be very near 1)
                error = np.sqrt((dx[xmin:xmax] * self.Spectrum.error[xmin:xmax]**2).sum() / dx[xmin:xmax].sum())
                return np.array([integ,error])
            else:
                return integ

            #OK = np.abs( fullmodel ) > threshold
            #integ = (self.spectofit[OK] * dx[OK]).sum()
            #error = np.sqrt((self.errspec[OK]**2 * dx[OK]).sum()/dx[OK].sum())
        else:
            if not hasattr(self.fitter,'integral'):
                raise AttributeError("The fitter %s does not have an integral implemented" % self.fittype)

            # the model considered here must NOT include the baseline!
            # if it does, you'll get the integral of the continuum
            #fullmodel = self.get_full_model(add_baseline=False)

            if self.Spectrum.xarr.cdelt() is not None:
                dx = np.median(dx)
                integ = self.fitter.integral(self.modelpars, dx=dx, **kwargs)
                if return_error:
                    if mycfg.WARN:
                        warn("WARNING: The computation of the error "
                             "on the integral is not obviously "
                             "correct or robust... it's just a guess.")
                    OK = self.model_mask(threshold=threshold, add_baseline=False)
                    error = np.sqrt((self.errspec[OK]**2).sum()) * dx
                    #raise NotImplementedError("We haven't written up correct error estimation for integrals of fits")
            else:
                integ = 0
                error = 0
                warn("An analytic integal could not be computed because the X-axis is irregular.  Try direct=True when integrating, or find a way to linearize the X-axis")
        if return_error:
            return integ,error
        else:
            return integ

    def model_mask(self, **kwargs):
        """
        Get a mask (boolean array) of the region where the fitted model is
        significant

        Parameters
        ----------
        threshold : 'auto' or 'error' or float
            The threshold to compare the model values to for selecting the mask
            region.

             * auto: uses `peak_fraction` times the model peak
             * error: use the spectrum error
             * float: any floating point number as an absolute threshold

        peak_fraction : float
            Parameter used if ``threshold=='auto'`` to determine fraction of
            model peak to set threshold at
        add_baseline : bool
            Add the fitted baseline to the model before comparing to threshold?

        Returns
        -------
        mask : `~numpy.ndarray`
            A boolean mask array with the same size as the spectrum, set to
            ``True`` where the fitted model has values above a specified
            threshold
        """
        return self._compare_to_threshold(**kwargs)

    def _compare_to_threshold(self, threshold='auto', peak_fraction=0.01,
                              add_baseline=False):
        """
        Identify pixels that are above some threshold
        """
        model = self.get_full_model(add_baseline=add_baseline)

        # auto-set threshold from some fraction of the model peak
        if threshold=='auto':
            threshold = peak_fraction * np.abs(model).max()
        elif threshold=='error':
            threshold = self.errspec

        OK = np.abs(model) > threshold

        return OK

    def get_model_xlimits(self, threshold='auto', peak_fraction=0.01,
                          add_baseline=False, unit='pixels', units=None):
        """
        Return the x positions of the first and last points at which the model
        is above some threshold

        Parameters
        ----------
        threshold : 'auto' or 'error' or float
            If 'auto', the threshold will be set to peak_fraction * the peak
            model value.
            If 'error', uses the error spectrum as the threshold
        peak_fraction : float
            ignored unless threshold == 'auto'
        add_baseline : bool
            Include the baseline when computing whether the model is above the
            threshold?  default FALSE.  Passed to get_full_model.
        units : str
            A valid unit type, e.g. 'pixels' or 'angstroms'
        """
        OK = self._compare_to_threshold(threshold=threshold,
                                        peak_fraction=peak_fraction,
                                        add_baseline=add_baseline)

        # find the first & last "True" values
        xpixmin = OK.argmax()
        xpixmax = len(OK) - OK[::-1].argmax() - 1

        if units is not None and unit =='pixels':
            # todo: deprecate
            unit = units

        if unit == 'pixels':
            return [xpixmin,xpixmax]
        else:
            return self.Spectrum.xarr[[xpixmin,xpixmax]].as_unit(units)

    def shift_pars(self, frame=None):
        """
        Shift the velocity / wavelength / frequency of the fitted parameters
        into a different frame

        Right now this only takes care of redshift and only if redshift is defined.
        It should be extended to do other things later
        """
        for ii,pi in enumerate(self.parinfo):
            for partype in ('shift','offset','velo'):
                if partype in str.lower(pi['parname']):
                    if frame is not None:
                        self.modelpars[ii] = self.Spectrum.xarr.x_in_frame(self.modelpars[ii], frame)

    def moments(self, fittype=None, **kwargs):
        """
        Return the moments

        see the :mod:`~pyspeckit.spectrum.moments` module

        Parameters
        ----------
        fittype : None or str
            The registered fit type to use for moment computation
        """
        if fittype is None:
            fittype = self.fittype
        return list(self.Registry.multifitters[fittype].moments(
            self.Spectrum.xarr[self.xmin:self.xmax],
            self.spectofit[self.xmin:self.xmax],  **kwargs))

    def button3action(self, event, debug=False, nwidths=1):
        """
        Disconnect the interactiveness
        Perform the fit (or die trying)
        Hide the guesses
        """
        if self.nclicks_b1 == 0:
            # there has been no selection
            # therefore, we assume *everything* is selected
            self.includemask[:] = True

        self.Spectrum.plotter.figure.canvas.mpl_disconnect(self.click)
        self.Spectrum.plotter.figure.canvas.mpl_disconnect(self.keyclick)
        if self.npeaks > 0:

            if hasattr(self, 'fitter'):
                self.guesses = self.fitter.parse_3par_guesses(self.guesses)
            else:
                # default fitter is a Gaussian, which has 3 parameters
                if len(self.guesses) % 3 != 0:
                    log.error("Default fitter is Gaussian, and there were "
                              "{0} guess parameters, which is not a "
                              "multiple of 3.".format(len(self.guesses)))

            log.info("{0} Guesses : {1}  X channel range: {2}-{3}"
                     .format(len(self.guesses), self.guesses, self.xmin,
                             self.xmax))

            self.multifit(use_window_limits=True)

            for p in self.button2plot + self.button1plot:
                p.set_visible(False)

        # disconnect interactive window (and more importantly, reconnect to
        # original interactive cmds)
        self.clear_all_connections()

    def copy(self, parent=None, registry=None):
        """
        Create a copy of the spectral fit - includes copies of the _full_model,
        the registry, the fitter, parinfo, modelpars, modelerrs, model, npeaks

        Parameters
        ----------
        parent : `pyspeckit.classes.Spectrum`
            A `~Spectrum` instance that is the parent of the specfit
            instance.  This needs to be specified at some point, but defaults
            to None to prevent overwriting a previous plot.
        """

        if registry is None:
            if hasattr(parent, 'Registry'):
                registry = parent.Registry
            else:
                # only make a copy if we're not already given a specific registry
                # to inherit
                copy.deepcopy(self.Registry)

        newspecfit = Specfit(parent, Registry=registry)
        newspecfit.parinfo = copy.deepcopy(self.parinfo)
        if newspecfit.parinfo is None:
            newspecfit.modelpars = None
            newspecfit.modelerrs = None
        else:
            newspecfit.modelpars = newspecfit.parinfo.values
            newspecfit.modelerrs = newspecfit.parinfo.errors
        newspecfit.includemask = self.includemask.copy()
        newspecfit.model = copy.copy(self.model)
        newspecfit.npeaks = self.npeaks
        if hasattr(self,'fitter'):
            newspecfit.fitter = copy.deepcopy(self.fitter)
            newspecfit.fitter.parinfo = newspecfit.parinfo
        if hasattr(self,'fullmodel'):
            newspecfit._full_model()

        # this is ridiculous, absurd, and infuriating...
        newspecfit.button2action = newspecfit.guesspeakwidth

        if parent is not None:
            newspecfit.Spectrum.plotter = parent.plotter
        else:
            newspecfit.Spectrum.plotter = None

        return newspecfit

    def __copy__(self):
        return self.copy(parent=self.Spectrum)

    def add_sliders(self, parlimitdict=None, **kwargs):
        """
        Add a Sliders window in a new figure appropriately titled

        Parameters
        ----------
        parlimitdict: dict
            Each parameter needs to have displayed limits; these are set in
            min-max pairs.  If this is left empty, the widget will try to guess
            at reasonable limits, but the guessing is not very sophisticated
            yet.

        .. todo:: Add a button in the navbar that makes this window pop up
        http://stackoverflow.com/questions/4740988/add-new-navigate-modes-in-matplotlib
        """

        if parlimitdict is None:
            # try to create a reasonable parlimit dict
            parlimitdict = {}

        for param in self.parinfo:
            if not param.parname in parlimitdict:
                if any( (x in param['parname'].lower() for x in ('shift','xoff')) ):
                    lower, upper = (self.Spectrum.xarr[self.includemask].min().value,
                                    self.Spectrum.xarr[self.includemask].max().value)
                elif any( (x in param['parname'].lower() for x in ('width','fwhm')) ):
                    xvalrange = (self.Spectrum.xarr[self.includemask].max().value -
                                 self.Spectrum.xarr[self.includemask].min().value)
                    lower,upper = (0,xvalrange)
                elif any( (x in param['parname'].lower() for x in ('amp','peak','height')) ):
                    datarange = self.spectofit.max() - self.spectofit.min()
                    lower,upper = (param['value']-datarange, param['value']+datarange)
                else:
                    lower = param['value'] * 0.1
                    upper = param['value'] * 10

                # override guesses with limits
                if param.limited[0]:
                    # nextafter -> next representable float
                    lower = np.nextafter(param.limits[0], param.limits[0]+1)
                if param.limited[1]:
                    upper = np.nextafter(param.limits[1], param.limits[1]-1)

                parlimitdict[param.parname] = (lower,upper)

        if hasattr(self,'fitter'):
            self.SliderWidget = widgets.FitterSliders(self,
                                                      self.Spectrum.plotter.figure,
                                                      npars=self.fitter.npars,
                                                      parlimitdict=parlimitdict,
                                                      **kwargs)
        else:
            log.error("Must have a fitter instantiated before creating sliders")

    def optimal_chi2(self, reduced=True, threshold='error', **kwargs):
        """
        Compute an "optimal" :math:`\\chi^2` statistic, i.e. one in which only pixels in
        which the model is statistically significant are included

        Parameters
        ----------
        reduced : bool
            Return the reduced :math:`\\chi^2`
        threshold : 'auto' or 'error' or float
            If 'auto', the threshold will be set to peak_fraction * the peak
            model value, where peak_fraction is a kwarg passed to
            get_model_xlimits reflecting the fraction of the model peak
            to consider significant
            If 'error', uses the error spectrum as the threshold
        kwargs : dict
            passed to :meth:`get_model_xlimits`

        Returns
        -------
        chi2 : float
            :math:`\\chi^2` statistic or reduced :math:`\\chi^2` statistic (:math:`\\chi^2/n`)

            .. math::
                   \\chi^2 = \\sum( (d_i - m_i)^2 / e_i^2 )
        """

        modelmask = self._compare_to_threshold(threshold=threshold, **kwargs)

        chi2 = np.sum((self.fullresiduals[modelmask]/self.errspec[modelmask])**2)

        if reduced:
            # vheight included here or not?  assuming it should be...
            return chi2/self.dof
        else:
            return chi2

    def get_pymc(self, **kwargs):
        """
        Create a pymc MCMC sampler from the current fitter.  Defaults to 'uninformative' priors

        `kwargs` are passed to the fitter's get_pymc method, with parameters defined below.

        Parameters
        ----------
        data : np.ndarray
        error : np.ndarray
        use_fitted_values : bool
            Each parameter with a measured error will have a prior defined by
            the Normal distribution with sigma = par.error and mean = par.value

        Examples
        --------

        >>> x = pyspeckit.units.SpectroscopicAxis(np.linspace(-10,10,50), unit='km/s')
        >>> e = np.random.randn(50)
        >>> d = np.exp(-np.asarray(x)**2/2.)*5 + e
        >>> sp = pyspeckit.Spectrum(data=d, xarr=x, error=np.ones(50)*e.std())
        >>> sp.specfit(fittype='gaussian')
        >>> MCuninformed = sp.specfit.get_pymc()
        >>> MCwithpriors = sp.specfit.get_pymc(use_fitted_values=True)
        >>> MCuninformed.sample(1000)
        >>> MCuninformed.stats()['AMPLITUDE0']
        >>> # WARNING: This will fail because width cannot be set <0, but it may randomly reach that...
        >>> # How do you define a likelihood distribution with a lower limit?!
        >>> MCwithpriors.sample(1000)
        >>> MCwithpriors.stats()['AMPLITUDE0']
        """
        if hasattr(self.fitter,'get_pymc'):
            return self.fitter.get_pymc(self.Spectrum.xarr, self.spectofit,
                                        self.errspec, **kwargs)
        else:
            raise AttributeError("Fitter %r does not have pymc implemented." % self.fitter)

    def get_emcee(self, nwalkers=None, **kwargs):
        """
        Get an emcee walker ensemble for the data & model using the current model type

        Parameters
        ----------
        data : np.ndarray
        error : np.ndarray
        nwalkers : int
            Number of walkers to use.  Defaults to 2 * self.fitters.npars

        Examples
        --------

        >>> import pyspeckit
        >>> x = pyspeckit.units.SpectroscopicAxis(np.linspace(-10,10,50), unit='km/s')
        >>> e = np.random.randn(50)
        >>> d = np.exp(-np.asarray(x)**2/2.)*5 + e
        >>> sp = pyspeckit.Spectrum(data=d, xarr=x, error=np.ones(50)*e.std())
        >>> sp.specfit(fittype='gaussian')
        >>> emcee_ensemble = sp.specfit.get_emcee()
        >>> p0 = emcee_ensemble.p0 * (np.random.randn(*emcee_ensemble.p0.shape) / 10. + 1.0)
        >>> pos,logprob,state = emcee_ensemble.run_mcmc(p0,100)
        """
        import emcee

        if hasattr(self.fitter,'get_emcee_ensemblesampler'):
            nwalkers = (self.fitter.npars * self.fitter.npeaks + self.fitter.vheight) * 2
            emc = self.fitter.get_emcee_ensemblesampler(self.Spectrum.xarr,
                                                        self.spectofit,
                                                        self.errspec, nwalkers)
            emc.nwalkers = nwalkers
            emc.p0 = np.array([self.parinfo.values] * emc.nwalkers)
            return emc

    def get_components(self, **kwargs):
        """
        If a model has been fitted, return the components of the model

        Parameters
        ----------
        kwargs are passed to self.fitter.components
        """
        if self.modelpars is not None:
            self.modelcomponents = self.fitter.components(self.Spectrum.xarr,
                    self.modelpars, **kwargs)

            return self.modelcomponents

    def measure_approximate_fwhm(self, threshold='error', emission=True,
                                 interpolate_factor=1, plot=False,
                                 grow_threshold=2, **kwargs):
        """
        Measure the FWHM of a fitted line

        This procedure is designed for multi-component *blended* lines; if the
        true FWHM is known (i.e., the line is well-represented by a single
        gauss/voigt/lorentz profile), use that instead.  Do not use this for
        multiple independently peaked profiles.

        This MUST be run AFTER a fit has been performed!

        Parameters
        ----------
        threshold : 'error' | float
            The threshold above which the spectrum will be interpreted as part
            of the line.  This threshold is applied to the *model*.  If it is
            'noise', self.error will be used.
        emission : bool
            Is the line absorption or emission?
        interpolate_factor : integer
            Magnification factor for determining sub-pixel FWHM.  If used,
            "zooms-in" by using linear interpolation within the line region
        plot : bool
            Overplot a line at the FWHM indicating the FWHM.  kwargs
            are passed to matplotlib.plot
        grow_threshold : int
            Minimum number of valid points.  If the total # of points above the
            threshold is <= to this number, it will be grown by 1 pixel on each side

        Returns
        -------
        The approximated FWHM, if it can be computed
        If there are <= 2 valid pixels, a fwhm cannot be computed
        """

        if threshold == 'error':
            threshold = self.Spectrum.error
            if np.all(self.Spectrum.error==0):
                threshold = 1e-3*self.Spectrum.data.max()

        if self.Spectrum.baseline.subtracted is False:
            data = self.Spectrum.data - self.Spectrum.baseline.basespec
        else:
            data = self.Spectrum.data * 1

        model = self.get_full_model(add_baseline=False)
        if np.count_nonzero(model) == 0:
            raise ValueError("The model is all zeros.  No FWHM can be "
                             "computed.")

        # can modify inplace because data is a copy of self.Spectrum.data
        if not emission:
            data *= -1
            model *= -1

        line_region = model > threshold
        if line_region.sum() == 0:
            raise ValueError("No valid data included in FWHM computation")
        if line_region.sum() <= grow_threshold:
            line_region[line_region.argmax()-1:line_region.argmax()+1] = True
            reverse_argmax = len(line_region) - line_region.argmax() - 1
            line_region[reverse_argmax-1:reverse_argmax+1] = True
            log.warning("Fewer than {0} pixels were identified as part of the fit."
                     " To enable statistical measurements, the range has been"
                     " expanded by 2 pixels including some regions below the"
                     " threshold.".format(grow_threshold))

        # determine peak (because data is neg if absorption, always use max)
        peak = data[line_region].max()

        xarr = self.Spectrum.xarr[line_region]
        xarr.make_dxarr()
        cd = xarr.dxarr.min()

        if interpolate_factor > 1:
            newxarr = units.SpectroscopicAxis(np.arange(xarr.min().value-cd.value,
                                                        xarr.max().value+cd.value,
                                                        cd.value /
                                                        float(interpolate_factor)
                                                       ),
                                              unit=xarr.unit,
                                              equivalencies=xarr.equivalencies
                                             )
            # load the metadata from xarr
            # newxarr._update_from(xarr)
            data = np.interp(newxarr,xarr,data[line_region])
            xarr = newxarr
        else:
            data = data[line_region]

        # need the peak location so we can find left/right half-max locations
        peakloc = data.argmax()

        hm_left = np.argmin(np.abs(data[:peakloc]-peak/2.))
        hm_right = np.argmin(np.abs(data[peakloc:]-peak/2.)) + peakloc

        deltax = xarr[hm_right]-xarr[hm_left]

        if plot:
            # for plotting, use a negative if absorption
            sign = 1 if emission else -1

            # shift with baseline if baseline is plotted
            if not self.Spectrum.baseline.subtracted:
                basespec = self.Spectrum.baseline.get_model(xarr)
                yoffleft = self.Spectrum.plotter.offset + basespec[hm_left]
                yoffright = self.Spectrum.plotter.offset + basespec[hm_right]
            else:
                yoffleft = yoffright = self.Spectrum.plotter.offset

            log.debug("peak={2} yoffleft={0} yoffright={1}".format(yoffleft, yoffright, peak))
            log.debug("hm_left={0} hm_right={1} xarr[hm_left]={2} xarr[hm_right]={3}".format(hm_left, hm_right, xarr[hm_left], xarr[hm_right]))

            self.Spectrum.plotter.axis.plot([xarr[hm_right].value,
                                             xarr[hm_left].value],
                                            np.array([sign*peak/2.+yoffleft,
                                                      sign*peak/2.+yoffright]),
                                            **kwargs)
            self.Spectrum.plotter.refresh()

        # debug print hm_left,hm_right,"FWHM: ",deltax
        # debug self.Spectrum.plotter.axis.plot(xarr,data,color='magenta')
        # debug self.Spectrum.plotter.refresh()
        # debug raise TheDead

        return deltax

    def _validate_parinfo(self, parinfo, mode='fix'):

        assert mode in ('fix','raise','check','guesses')

        any_out_of_range = []

        for param in parinfo:
            if (param.limited[0] and (param.value < param.limits[0])):
                if (np.allclose(param.value, param.limits[0])):
                    # nextafter -> next representable float
                    if mode in ('fix', 'guesses'):
                        warn("{0} is less than the lower limit {1}, but very close."
                             " Converting to {1}+ULP".format(param.value,
                                                             param.limits[0]))
                        param.value = np.nextafter(param.limits[0], param.limits[0]+1)
                    elif mode == 'raise':
                        raise ValueError("{0} is less than the lower limit {1}, but very close."
                                         .format(param.value, param.limits[1]))
                    elif mode == 'check':
                        any_out_of_range.append("lt:close",)
                else:
                    raise ValueError("{0} is less than the lower limit {1}"
                                     .format(param.value, param.limits[0]))
            elif mode == 'check':
                any_out_of_range.append(False)

            if (param.limited[1] and (param.value > param.limits[1])):
                if (np.allclose(param.value, param.limits[1])):
                    if mode in ('fix', 'guesses'):
                        param.value = np.nextafter(param.limits[1], param.limits[1]-1)
                        warn("{0} is greater than the upper limit {1}, but very close."
                             " Converting to {1}-ULP".format(param.value,
                                                             param.limits[1]))
                    elif mode == 'raise':
                        raise ValueError("{0} is greater than the upper limit {1}, but very close."
                                         .format(param.value, param.limits[1]))
                    elif mode == 'check':
                        any_out_of_range.append("gt:close")
                else:
                    raise ValueError("{0} is greater than the upper limit {1}"
                                     .format(param.value, param.limits[1]))
            elif mode == 'check':
                any_out_of_range.append(False)

        if mode == 'guesses':
            return parinfo.values
        return any_out_of_range
