"""
=============================
Generic SpectralModel wrapper
=============================
.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>

Module API
^^^^^^^^^^
"""
import numpy as np
from pyspeckit.mpfit import mpfit,mpfitException
from pyspeckit.spectrum.parinfo import ParinfoList,Parinfo
import copy
from astropy import log
import matplotlib.cbook as mpcb
from . import fitter
from . import mpfit_messages
from pyspeckit.specwarnings import warn
from pyspeckit.spectrum.units import SpectroscopicAxis
import itertools
import operator
import six
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict
except ImportError:
    warn("OrderedDict is required for modeling.  "
         "If you have python <2.7, install the ordereddict module.")

# define the allowed guess types and the order in which they are received
valid_guess_types = ('amplitude', 'center', 'width')

class SpectralModel(fitter.SimpleFitter):
    """
    A wrapper class for a spectra model.  Includes internal functions to
    generate multi-component models, annotations, integrals, and individual
    components.  The declaration can be complex, since you should name
    individual variables, set limits on them, set the units the fit will be
    performed in, and set the annotations to be used.  Check out some
    of the hyperfine codes (hcn, n2hp) for examples.
    """

    def __init__(self, modelfunc, npars,
                 shortvarnames=("A","\\Delta x","\\sigma"),
                 fitunit=None,
                 centroid_par=None,
                 fwhm_func=None,
                 fwhm_pars=None,
                 integral_func=None,
                 use_lmfit=False,
                 guess_types=('amplitude', 'center', 'width'),
                 **kwargs):
        """
        Spectral Model Initialization

        Create a Spectral Model class for data fitting

        Parameters
        ----------
        modelfunc : function
            the model function to be fitted.  Should take an X-axis
            (spectroscopic axis) as an input followed by input parameters.
            Returns an array with the same shape as the input X-axis
        npars : int
            number of parameters required by the model
        use_lmfit: bool
            Use lmfit instead of mpfit to do the fitting
        parnames : list (optional)
            a list or tuple of the parameter names
        parvalues : list (optional)
            the initial guesses for the input parameters (defaults to ZEROS)
        parlimits :  list (optional)
            the upper/lower limits for each variable (defaults to ZEROS)
        parfixed  : list (optional)
            Can declare any variables to be fixed (defaults to ZEROS)
        parerror  : list (optional)
            technically an output parameter.  Specifying it here will have no
            effect. (defaults to ZEROS)
        partied   : list (optional)
            not the past tense of party.  Can declare, via text, that
            some parameters are tied to each other.  Defaults to zeros like the
            others, but it's not clear if that's a sensible default
        fitunit : str (optional)
            convert X-axis to these units before passing to model
        parsteps : list (optional)
            minimum step size for each paremeter (defaults to ZEROS)
        npeaks   : list (optional)
            default number of peaks to assume when fitting (can be overridden)
        shortvarnames : list (optional)
            TeX names of the variables to use when annotating
        amplitude_types : tuple
            A tuple listing the types of the different parameters when guessing.
            The valid values are 'amplitude', 'width', and 'center'.  These are
            handled by parse_3par_guesses, which translate these into input
            guess lists for the fitter.  For a "standard" 3-parameter Gaussian
            fitter, nothing changes, but for other models that have more than
            3 parameters, some translation is needed.

        Returns
        -------
        A tuple containing (model best-fit parameters, the model, parameter
        errors, chi^2 value)
        """

        self.modelfunc = modelfunc
        if self.__doc__ is None:
            self.__doc__ = modelfunc.__doc__
        elif modelfunc.__doc__ is not None:
            self.__doc__ += modelfunc.__doc__
        self.npars = npars
        self.default_npars = npars
        self.fitunit = fitunit

        # this needs to be set once only
        self.shortvarnames = shortvarnames

        self.default_parinfo = None
        self.default_parinfo, kwargs = self._make_parinfo(**kwargs)
        self.parinfo = copy.copy(self.default_parinfo)

        self.modelfunc_kwargs = kwargs

        self.use_lmfit = use_lmfit

        # default name of parameter that represents the profile centroid
        self.centroid_par = centroid_par

        # FWHM function and parameters
        self.fwhm_func = fwhm_func
        self.fwhm_pars = fwhm_pars

        # analytic integral function
        self.integral_func = integral_func

        for gt in guess_types:
            if not isinstance(gt, float) and not any(g in gt for g in valid_guess_types):
                raise ValueError("Guess type must be one of {0} or a float"
                                 .format(valid_guess_types))
        self.guess_types = guess_types

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

    def __call__(self, *args, **kwargs):

        log.debug("Fitter called with args={0} and kwargs={1}".format(args, kwargs))
        use_lmfit = kwargs.pop('use_lmfit') if 'use_lmfit' in kwargs else self.use_lmfit
        if use_lmfit:
            return self.lmfitter(*args,**kwargs)
        return self.fitter(*args,**kwargs)

    @property
    def npeaks(self):
        return int(self._npeaks)

    @npeaks.setter
    def npeaks(self, value):
        if int(value) != value:
            raise ValueError("npeaks must be an integer")
        self._npeaks = int(value)

    def make_parinfo(self, **kwargs):
        return self._make_parinfo(**kwargs)[0]

    def _make_parinfo(self, params=None, parnames=None, parvalues=None,
                      parlimits=None, parlimited=None, parfixed=None,
                      parerror=None, partied=None, fitunit=None,
                      parsteps=None, npeaks=1, parinfo=None, names=None,
                      values=None, limits=None, limited=None, fixed=None,
                      error=None, tied=None, steps=None, negamp=None,
                      limitedmin=None, limitedmax=None, minpars=None,
                      maxpars=None, vheight=False, debug=False, **kwargs):
        """
        Generate a `ParinfoList` that matches the inputs

        This code is complicated - it can take inputs in a variety of different
        forms with different priority.  It will return a `ParinfoList` (and
        therefore must have values within parameter ranges)

        """
        log.debug("BEGIN _make_parinfo")

        # for backwards compatibility - partied = tied, etc.
        locals_dict = locals()
        for varname in str.split("parnames,parvalues,parsteps,parlimits,parlimited,parfixed,parerror,partied",","):
            shortvarname = varname.replace("par","")
            if locals_dict.get(shortvarname) is not None and locals_dict.get(varname) is not None:
                raise ValueError("Cannot specify both {0} and {1}".format(varname, shortvarname))

        input_pardict = {k: locals_dict.get(k)
                         for k in str.split("parnames,parvalues,parsteps,parlimits,parlimited,parfixed,parerror,partied",",")}
        _tip = {'par'+k: locals_dict.get(k)
                for k in str.split("names,values,steps,limits,limited,fixed,error,tied",",")
                if locals_dict.get(k)
               }
        input_pardict.update(_tip)

        if params is not None and parvalues is not None:
            raise ValueError("parvalues and params both specified; they're redundant so that's not allowed.")
        elif params is not None and parvalues is None:
            input_pardict['parvalues'] = params
        log.debug("Parvalues = {0}, npeaks = {1}".format(input_pardict['parvalues'], npeaks))

        # this is used too many damned times to keep referencing a dict.
        parnames = input_pardict['parnames']
        parlimited = input_pardict['parlimited']
        parlimits = input_pardict['parlimits']
        parvalues = input_pardict['parvalues']

        if parnames is not None:
            self.parnames = parnames
        elif parnames is None and hasattr(self,'parnames') and self.parnames is not None:
            parnames = self.parnames
        elif self.default_parinfo is not None and parnames is None:
            parnames = [p['parname'] for p in self.default_parinfo]

        input_pardict['parnames'] = parnames

        assert input_pardict['parnames'] is not None

        if limitedmin is not None:
            if limitedmax is not None:
                parlimited = list(zip(limitedmin,limitedmax))
            else:
                parlimited = list(zip(limitedmin,(False,)*len(parnames)))
        elif limitedmax is not None:
            parlimited = list(zip((False,)*len(parnames),limitedmax))
        elif self.default_parinfo is not None and parlimited is None:
            parlimited = [p['limited'] for p in self.default_parinfo]

        input_pardict['parlimited'] = parlimited

        if minpars is not None:
            if maxpars is not None:
                parlimits = list(zip(minpars,maxpars))
            else:
                parlimits = list(zip(minpars,(False,)*len(parnames)))
        elif maxpars is not None:
            parlimits = list(zip((False,)*len(parnames),maxpars))
        elif limits is not None:
            parlimits = limits
        elif self.default_parinfo is not None and parlimits is None:
            parlimits = [p['limits'] for p in self.default_parinfo]

        input_pardict['parlimits'] = parlimits

        self.npeaks = int(npeaks)

        # the height / parvalue popping needs to be done before the temp_pardict is set in order to make sure
        # that the height guess isn't assigned to the amplitude
        self.vheight = vheight
        if ((vheight and len(self.parinfo) == self.default_npars and
             len(parvalues) == self.default_npars + 1)):
            # if the right number of parameters are passed, the first is the height
            self.parinfo = [{'n':0, 'value':parvalues.pop(0), 'limits':(0,0),
                             'limited': (False,False), 'fixed':False, 'parname':'HEIGHT',
                             'error': 0, 'tied':""}]
        elif vheight and len(self.parinfo) == self.default_npars and len(parvalues) == self.default_npars:
            # if you're one par short, guess zero
            self.parinfo = [{
                'n':0, 'value': 0, 'limits':(0,0), 'limited': (False,False),
                'fixed':False, 'parname':'HEIGHT', 'error': 0, 'tied':""
            }]
        elif vheight and len(self.parinfo) == self.default_npars+1 and len(parvalues) == self.default_npars+1:
            # the right numbers are passed *AND* there is already a height param
            self.parinfo = [{
                'n':0, 'value':parvalues.pop(0), 'limits':(0,0),
                'limited': (False,False), 'fixed': False, 'parname':'HEIGHT',
                'error': 0, 'tied':""
            }]
            #heightparnum = (i for i,s in self.parinfo if 'HEIGHT' in s['parname'])
            #for hpn in heightparnum:
            #    self.parinfo[hpn]['value'] = parvalues[0]
        elif vheight:
            raise ValueError('VHEIGHT is specified but a case was found that did not allow it to be included.')
        else:
            self.parinfo = []

        log.debug("After VHEIGHT parse len(parinfo): %i   vheight: %s" % (len(self.parinfo), vheight))


        # this is a clever way to turn the parameter lists into a dict of lists
        # clever = hard to read
        temp_pardict = OrderedDict([(varname, np.zeros(self.npars*self.npeaks,
                                                       dtype='bool'))
                                    if input_pardict.get(varname) is None else
                                    (varname, list(input_pardict.get(varname)))
            for varname in str.split("parnames,parvalues,parsteps,parlimits,parlimited,parfixed,parerror,partied",",")])
        temp_pardict['parlimits'] = parlimits if parlimits is not None else [(0,0)] * (self.npars*self.npeaks)
        temp_pardict['parlimited'] = parlimited if parlimited is not None else [(False,False)] * (self.npars*self.npeaks)
        for k,v in temp_pardict.items():
            if (self.npars*self.npeaks) / len(v) > 1:
                n_components = ((self.npars*self.npeaks) / len(v))
                if n_components != int(n_components):
                    raise ValueError("The number of parameter values is not a "
                                     "multiple of the number of allowed "
                                     "parameters.")
                temp_pardict[k] = list(v) * int(n_components)

        # generate the parinfo dict
        # note that 'tied' must be a blank string (i.e. ""), not False, if it is not set
        # parlimited, parfixed, and parlimits are all two-element items (tuples or lists)
        self.parinfo += [{'n':ii+self.npars*jj+vheight,
                          'value':float(temp_pardict['parvalues'][ii+self.npars*jj]),
                          'step':temp_pardict['parsteps'][ii+self.npars*jj],
                          'limits':temp_pardict['parlimits'][ii+self.npars*jj],
                          'limited':temp_pardict['parlimited'][ii+self.npars*jj],
                          'fixed':temp_pardict['parfixed'][ii+self.npars*jj],
                          'parname':temp_pardict['parnames'][ii].upper()+"%0i" % int(jj),
                          'error':float(temp_pardict['parerror'][ii+self.npars*jj]),
                          'tied':temp_pardict['partied'][ii+self.npars*jj] if temp_pardict['partied'][ii+self.npars*jj] else ""}
            for jj in range(self.npeaks)
            for ii in range(self.npars) ] # order matters!

        log.debug("After Generation step len(parinfo): %i   vheight: %s "
                  "parinfo: %s" % (len(self.parinfo), vheight, self.parinfo))

        if debug > True:
            import pdb; pdb.set_trace()

        # special keyword to specify emission/absorption lines
        if negamp is not None:
            if negamp:
                for p in self.parinfo:
                    if 'AMP' in p['parname']:
                        p['limited'] = (p['limited'][0], True)
                        p['limits']  = (p['limits'][0],  0)
            else:
                for p in self.parinfo:
                    if 'AMP' in p['parname']:
                        p['limited'] = (True, p['limited'][1])
                        p['limits']  = (0, p['limits'][1])

        # This is effectively an override of all that junk above (3/11/2012)
        # Much of it is probably unnecessary, but it was easier to do this than
        # rewrite the above
        self.parinfo = ParinfoList([Parinfo(p) for p in self.parinfo])

        # New feature: scaleability
        for par in self.parinfo:
            if par.parname.lower().strip('0123456789') in ('amplitude','amp'):
                par.scaleable = True

        log.debug("Parinfo has been set: {0}".format(self.parinfo))
        log.debug("kwargs {0} were passed.".format(kwargs))

        assert self.parinfo != []

        return self.parinfo, kwargs

    def n_modelfunc(self, pars=None, debug=False, **kwargs):
        """
        Simple wrapper to deal with N independent peaks for a given spectral model
        """
        if pars is None:
            pars = self.parinfo
        elif not isinstance(pars, ParinfoList):
            try:
                partemp = copy.copy(self.parinfo)
                partemp._from_Parameters(pars)
                pars = partemp
            except AttributeError:
                log.log(5, "Reading pars {0} as LMPar failed.".format(pars))
                if debug > 1:
                    import pdb; pdb.set_trace()
        if hasattr(pars,'values'):
            # important to treat as Dictionary, since lmfit params & parinfo both have .items
            parnames,parvals = list(zip(*list(pars.items())))
            parnames = [p.lower() for p in parnames]
            parvals = [p.value for p in parvals]
        else:
            parvals = list(pars)

        if np.any(np.isnan(parvals)):
            raise ValueError("A parameter is NaN.  Unless you gave a NaN "
                             "value directly, this is a bug and should be "
                             "reported.  If you specified a NaN parameter, "
                             "don't do that.")

        log.debug("pars to n_modelfunc: {0}, parvals:{1}".format(pars, parvals))
        def L(x):
            if hasattr(x, 'value') and not hasattr(x, 'x_to_coord'):
                x = SpectroscopicAxis(x)
            v = np.zeros(len(x))
            if self.vheight:
                v += parvals[0]
            # use len(pars) instead of self.npeaks because we want this to work
            # independent of the current best fit
            for jj in range(int((len(parvals)-self.vheight)/self.npars)):
                lower_parind = jj*self.npars+self.vheight
                upper_parind = (jj+1)*self.npars+self.vheight
                v += self.modelfunc(x, *parvals[lower_parind:upper_parind], **kwargs)

            return v
        return L

    def mpfitfun(self,x,y,err=None):
        """
        Wrapper function to compute the fit residuals in an mpfit-friendly format
        """
        if err is None:
            def f(p,fjac=None):
                residuals = (y-self.n_modelfunc(p, **self.modelfunc_kwargs)(x))
                return [0,residuals]
        else:
            def f(p,fjac=None):
                residuals = (y-self.n_modelfunc(p, **self.modelfunc_kwargs)(x))/err
                return [0,residuals]
        return f

    def lmfitfun(self,x,y,err=None,debug=False):
        """
        Wrapper function to compute the fit residuals in an lmfit-friendly format
        """
        def f(p):
            #pars = [par.value for par in p.values()]
            kwargs = {}
            kwargs.update(self.modelfunc_kwargs)

            log.debug("Pars, kwarg keys: {0},{1}".format(p,list(kwargs.keys())))
            if err is None:
                return (y-self.n_modelfunc(p,**kwargs)(x))
            else:
                return (y-self.n_modelfunc(p,**kwargs)(x))/err
        return f

    def lmfitter(self, xax, data, err=None, parinfo=None, quiet=True, debug=False, **kwargs):
        """
        Use lmfit instead of mpfit to do the fitting

        Parameters
        ----------
        xax : SpectroscopicAxis
            The X-axis of the spectrum
        data : ndarray
            The data to fit
        err : ndarray (optional)
            The error on the data.  If unspecified, will be uniform unity
        parinfo : ParinfoList
            The guesses, parameter limits, etc.  See
            `pyspeckit.spectrum.parinfo` for details
        quiet : bool
            If false, print out some messages about the fitting

        """
        try:
            import lmfit
        except ImportError as e:
            raise ImportError("Could not import lmfit, try using mpfit instead.")

        log.debug("lmfit called with parinfo=\n{0}".format(parinfo))

        self.xax = xax # the 'stored' xax is just a link to the original
        if hasattr(xax,'convert_to_unit') and self.fitunit is not None:
            # some models will depend on the input units.  For these, pass in an X-axis in those units
            # (gaussian, voigt, lorentz profiles should not depend on units.  Ammonia, formaldehyde,
            # H-alpha, etc. should)
            xax = copy.copy(xax)
            xax.convert_to_unit(self.fitunit, quiet=quiet)
        elif self.fitunit is not None:
            raise TypeError("X axis does not have a convert method")

        if np.any(np.isnan(data)) or np.any(np.isinf(data)):
            err[np.isnan(data) + np.isinf(data)] = np.inf
            data[np.isnan(data) + np.isinf(data)] = 0
        if np.any(np.isnan(err)):
            raise ValueError("One or more of the error values is NaN."
                             "  This is not allowed.  Errors can be infinite "
                             "(which is equivalent to giving zero weight to "
                             "a data point), but otherwise they must be positive "
                             "floats.")
        elif np.any(err<0):
            raise ValueError("At least one error value is negative, which is "
                             "not allowed as negative errors are not "
                             "meaningful in the optimization process.")


        if parinfo is None:
            parinfo, kwargs = self._make_parinfo(debug=debug, **kwargs)
            log.debug("Parinfo created from _make_parinfo: {0}".format(parinfo))

        LMParams = parinfo.as_Parameters()
        log.debug("LMParams: "+"\n".join([repr(p) for p in list(LMParams.values())]))
        log.debug("parinfo:  {0}".format(parinfo))
        log.debug("BEGIN MINIMIZER")
        minimizer = lmfit.minimize(self.lmfitfun(xax,np.array(data),err,debug=debug),LMParams,**kwargs)
        log.debug("END MINIMIZER")
        if not quiet:
            log.info("There were %i function evaluations" % (minimizer.nfev))
        #modelpars = [p.value for p in parinfo.values()]
        #modelerrs = [p.stderr for p in parinfo.values() if p.stderr is not None else 0]

        self.LMParams = minimizer.params
        # Force consistency w/earlier versions of lmfit: if error == 0 exactly,
        # change it to None
        for par in self.LMParams:
            if hasattr(par, 'stderr') and par.stderr == 0:
                #assert minimizer.ier == 4
                par.stderr = None

        self.parinfo._from_Parameters(self.LMParams)
        log.debug("LMParams: {0}".format(self.LMParams))
        log.debug("parinfo: {0}".format(parinfo))


        self.mp = minimizer
        self.mpp = self.parinfo.values
        self.mpperr = self.parinfo.errors
        self.mppnames = self.parinfo.names
        modelkwargs = {}
        modelkwargs.update(self.modelfunc_kwargs)
        self.model = self.n_modelfunc(self.parinfo, **modelkwargs)(xax)
        if hasattr(minimizer,'chisqr'):
            chi2 = minimizer.chisqr
        else:
            try:
                chi2 = (((data-self.model)/err)**2).sum()
            except TypeError:
                chi2 = ((data-self.model)**2).sum()
        if np.isnan(chi2):
            warn("Warning: chi^2 is nan")

        if hasattr(self.mp,'ier') and self.mp.ier not in [1,2,3,4]:
            log.warning("Fitter failed: %s, %s" % (self.mp.message, self.mp.lmdif_message))

        return self.mpp,self.model,self.mpperr,chi2

    def fitter(self, xax, data, err=None, quiet=True, veryverbose=False,
               debug=False, parinfo=None, **kwargs):
        """
        Run the fitter using mpfit.

        kwargs will be passed to _make_parinfo and mpfit.

        Parameters
        ----------
        xax : SpectroscopicAxis
            The X-axis of the spectrum
        data : ndarray
            The data to fit
        err : ndarray (optional)
            The error on the data.  If unspecified, will be uniform unity
        parinfo : ParinfoList
            The guesses, parameter limits, etc.  See
            `pyspeckit.spectrum.parinfo` for details
        quiet : bool
            pass to mpfit.  If False, will print out the parameter values for
            each iteration of the fitter
        veryverbose : bool
            print out a variety of mpfit output parameters
        debug : bool
            raise an exception (rather than a warning) if chi^2 is nan
        """

        if parinfo is None:
            parinfo, kwargs = self._make_parinfo(debug=debug, **kwargs)
        else:
            log.debug("Using user-specified parinfo dict")
            # clean out disallowed kwargs (don't want to pass them to mpfit)
            #throwaway, kwargs = self._make_parinfo(debug=debug, **kwargs)

        self.xax = xax # the 'stored' xax is just a link to the original
        if hasattr(xax,'as_unit') and self.fitunit is not None:
            # some models will depend on the input units.  For these, pass in an X-axis in those units
            # (gaussian, voigt, lorentz profiles should not depend on units.  Ammonia, formaldehyde,
            # H-alpha, etc. should)
            xax = copy.copy(xax)
            # xax.convert_to_unit(self.fitunit, quiet=quiet)
            xax = xax.as_unit(self.fitunit, quiet=quiet, **kwargs)
        elif self.fitunit is not None:
            raise TypeError("X axis does not have a convert method")

        if np.any(np.isnan(data)) or np.any(np.isinf(data)):
            err[np.isnan(data) + np.isinf(data)] = np.inf
            data[np.isnan(data) + np.isinf(data)] = 0

        if np.any(np.isnan(err)):
            raise ValueError("One or more of the error values is NaN."
                             "  This is not allowed.  Errors can be infinite "
                             "(which is equivalent to giving zero weight to "
                             "a data point), but otherwise they must be positive "
                             "floats.")
        elif np.any(err<0):
            raise ValueError("At least one error value is negative, which is "
                             "not allowed as negative errors are not "
                             "meaningful in the optimization process.")

        for p in parinfo: log.debug( p )
        log.debug( "\n".join(["%s %i: tied: %s value: %s" % (p['parname'],p['n'],p['tied'],p['value']) for p in parinfo]) )

        mp = mpfit(self.mpfitfun(xax,data,err),parinfo=parinfo,quiet=quiet,debug=debug,**kwargs)
        mpp = mp.params
        if mp.perror is not None:
            mpperr = mp.perror
        else:
            mpperr = mpp*0
        chi2 = mp.fnorm

        if mp.status == 0:
            if "parameters are not within PARINFO limits" in mp.errmsg:
                log.warning(parinfo)
            raise mpfitException(mp.errmsg)

        for i,(p,e) in enumerate(zip(mpp,mpperr)):
            self.parinfo[i]['value'] = p
            # for consistency w/lmfit, and because it makes more sense, errors
            # of 0 will instead be None
            self.parinfo[i]['error'] = e if (e != 0 or mp.status != 4) else None

        # sanity check: if status==4, errors could not be computed
        # Apparently some parameters can have errors estimated even if all can't?
        #if mp.status == 4:
        #    assert all([self.parinfo[ii]['error'] is None
        #                for ii in range(len(mpp))])

        if veryverbose:
            log.info("Fit status: {0}".format(mp.status))
            log.info("Fit error message: {0}".format(mp.errmsg))
            log.info("Fit message: {0}".format(mpfit_messages[mp.status]))
            for i,p in enumerate(mpp):
                log.info("{0}: {1} +/- {2}".format(self.parinfo[i]['parname'],
                                                    p,mpperr[i]))
            log.info("Chi2: {0} Reduced Chi2: {1}  DOF:{2}".format(mp.fnorm,
                                                                   mp.fnorm/(len(data)-len(mpp)),
                                                                   len(data)-len(mpp)))

        self.mp = mp
        self.mpp = self.parinfo.values
        self.mpperr = self.parinfo.errors
        self.mppnames = self.parinfo.names
        self.model = self.n_modelfunc(self.parinfo,**self.modelfunc_kwargs)(xax)
        log.debug("Modelpars: {0}".format(self.mpp))
        if np.isnan(chi2):
            if debug:
                raise ValueError("Error: chi^2 is nan")
            else:
                log.warning("Warning: chi^2 is nan")
        return mpp,self.model,mpperr,chi2

    def slope(self, xinp):
        """
        Find the local slope of the model at location x
        (x must be in xax's units)
        """
        if hasattr(self, 'model'):
            dm = np.diff(self.model)
            # convert requested x to pixels
            xpix = self.xax.x_to_pix(xinp)
            dmx = np.average(dm[xpix-1:xpix+1])
            if np.isfinite(dmx):
                return dmx
            else:
                return 0

    def annotations(self, shortvarnames=None, debug=False):
        """
        Return a list of TeX-formatted labels

        The values and errors are formatted so that only the significant digits
        are displayed.  Rounding is performed using the decimal package.

        Parameters
        ----------
        shortvarnames : list
            A list of variable names (tex is allowed) to include in the
            annotations.  Defaults to self.shortvarnames

        Examples
        --------
        >>> # Annotate a Gaussian
        >>> sp.specfit.annotate(shortvarnames=['A','\\Delta x','\\sigma'])
        """
        from decimal import Decimal # for formatting
        svn = self.shortvarnames if shortvarnames is None else shortvarnames
        # if pars need to be replicated....
        if len(svn) < self.npeaks*self.npars:
            svn = svn * self.npeaks

        parvals = self.parinfo.values
        parerrs = self.parinfo.errors

        loop_list = [(parvals[ii+jj*self.npars+self.vheight],
                      parerrs[ii+jj*self.npars+self.vheight],
                      svn[ii+jj*self.npars],
                      self.parinfo.fixed[ii+jj*self.npars+self.vheight],
                      jj)
                     for jj in range(self.npeaks) for ii in range(self.npars)]

        label_list = []
        for (value, error, varname, fixed, varnumber) in loop_list:
            log.debug(", ".join([str(x) for x in (value, error, varname, fixed, varnumber)]))
            if None in (value, error):
                label = "{0}({1})=None".format(varname, varnumber)
            elif fixed or error==0:
                label = ("$%s(%i)$=%8s" % (varname, varnumber,
                        Decimal("%g" % value).quantize(Decimal("%0.6g" % (value)))))
            else:
                label = ("$%s(%i)$=%8s $\\pm$ %8s" % (varname, varnumber,
                    Decimal("%g" % value).quantize(Decimal("%0.2g" % (min(np.abs([value,error]))))),
                    Decimal("%g" % error).quantize(Decimal("%0.2g" % (error))),))
            label_list.append(label)

        labels = tuple(mpcb.flatten(label_list))
        return labels

    def components(self, xarr, pars, **kwargs):
        """
        Return a numpy ndarray of shape [npeaks x modelshape] of the
        independent components of the fits
        """

        modelcomponents = np.array(
            [self.modelfunc(xarr,
                *pars[i*self.npars:(i+1)*self.npars],
                **dict(list(self.modelfunc_kwargs.items())+list(kwargs.items())))
            for i in range(self.npeaks)])

        if len(modelcomponents.shape) == 3:
            newshape = [modelcomponents.shape[0]*modelcomponents.shape[1], modelcomponents.shape[2]]
            modelcomponents = np.reshape(modelcomponents, newshape)

        return modelcomponents

    def integral(self, modelpars, dx=None, **kwargs):
        """
        Extremely simple integrator:
        IGNORES modelpars;
        just sums self.model
        """
        if not hasattr(self,'model'):
            raise ValueError("Must fit (or compute) a model before computing"
                             " its integral.")

        if dx is not None:
            return (self.model*dx).sum()
        else:
            return self.model.sum()

    def analytic_integral(self, modelpars=None, npeaks=None, npars=None):
        """
        Placeholder for analyic integrals; these must be defined for individual models
        """
        if self.integral_func is None:
            raise NotImplementedError("Analytic integrals must be implemented independently for each model type")

        # all of these parameters are allowed to be overwritten
        if modelpars is None:
            modelpars = self.parinfo.values
        if npeaks is None:
            npeaks = self.npeaks
        if npars is None:
            npars = self.npars

        return np.sum([
            self.integral_func(modelpars[npars*ii:npars*(1+ii)])
            for ii in range(npeaks)])

    def component_integrals(self, xarr, dx=None):
        """
        Compute the integrals of each component
        """
        components = self.components(xarr, self.parinfo.values)

        if dx is None:
            dx = 1

        integrals = [com.sum()*dx for com in components]
        return integrals

    def analytic_fwhm(self, parinfo=None):
        """
        Return the FWHMa of the model components *if* a fwhm_func has been
        defined

        Done with incomprehensible list comprehensions instead of nested for
        loops... readability sacrificed for speed and simplicity.  This is
        unpythonic.
        """
        if self.fwhm_func is None and self.fwhm_pars is None:
            raise TypeError("fwhm_func not implemented for model %s" % self.__name__)

        if parinfo is None:
            parinfo = self.parinfo

        fwhm = [self.fwhm_func(
                *[self.parinfo[str.upper(p+'%i' % n)] for p in self.fwhm_pars]
                )
                for n in range(self.npeaks)]
        return fwhm

    def analytic_centroids(self, centroidpar=None):
        """
        Return the *analytic* centroids of the model components

        Parameters
        ----------
        centroidpar : None or string
            The name of the parameter in the fit that represents the centroid
            *some models have default centroid parameters - these will be used
            if centroidpar is unspecified*

        Returns
        -------
        List of the centroid values (even if there's only 1)
        """
        if centroidpar is None:
            centroidpar = self.centroid_par

        centr = [par.value
                for par in self.parinfo
                if str.upper(centroidpar) in par.parname]

        return centr


    def computed_centroid(self, xarr=None):
        """
        Return the *computed* centroid of the model

        Parameters
        ----------
        xarr : None or np.ndarray
            The X coordinates of the model over which the centroid should be
            computed.  If unspecified, the centroid will be in pixel units
        """
        if not hasattr(self, 'model'):
            raise ValueError("Must fit (or compute) a model before measuring "
                             "its centroid")
        if xarr is None:
            xarr = np.arange(self.model.size)

        centr = (self.model*xarr).sum() / self.model.sum()
        return centr

    def logp(self, xarr, data, error, pars=None):
        """
        Return the log probability of the model.  If the parameter is out of
        range, return -inf
        """
        if pars is None:
            pars = self.parinfo
        else:
            parinfo = copy.copy(self.parinfo)
            for value,parameter in zip(pars,parinfo):
                try:
                    parameter.value = value
                except ValueError:
                    return -np.inf
        model = self.n_modelfunc(pars, **self.modelfunc_kwargs)(xarr)

        difference = np.abs(data-model)

        # prob = 1/(2*np.pi)**0.5/error * exp(-difference**2/(2.*error**2))

        #logprob = np.log(1./(2.*np.pi)**0.5/error) * (-difference**2/(2.*error**2))
        logprob = (-difference**2/(2.*error**2))

        totallogprob = np.sum(logprob)

        return totallogprob

    def get_emcee_sampler(self, xarr, data, error, **kwargs):
        """
        Get an emcee walker for the data & model

        Parameters
        ----------
        xarr : pyspeckit.units.SpectroscopicAxis
        data : np.ndarray
        error : np.ndarray

        Examples
        --------

        >>> import pyspeckit
        >>> x = pyspeckit.units.SpectroscopicAxis(np.linspace(-10,10,50), unit='km/s')
        >>> e = np.random.randn(50)
        >>> d = np.exp(-np.asarray(x)**2/2.)*5 + e
        >>> sp = pyspeckit.Spectrum(data=d, xarr=x, error=np.ones(50)*e.std())
        >>> sp.specfit(fittype='gaussian')
        >>> emcee_sampler = sp.specfit.fitter.get_emcee_sampler(sp.xarr, sp.data, sp.error)
        >>> p0 = sp.specfit.parinfo
        >>> emcee_sampler.run_mcmc(p0,100)
        """
        try:
            import emcee
        except ImportError:
            return

        def probfunc(pars):
            return self.logp(xarr, data, error, pars=pars)

        raise NotImplementedError("emcee's metropolis-hastings sampler is not implemented; use pymc")
        sampler = emcee.MHSampler(self.npars*self.npeaks+self.vheight, probfunc, **kwargs)

        return sampler

    def get_emcee_ensemblesampler(self, xarr, data, error, nwalkers, **kwargs):
        """
        Get an emcee walker ensemble for the data & model

        Parameters
        ----------
        data : np.ndarray
        error : np.ndarray
        nwalkers : int
            Number of walkers to use

        Examples
        --------

        >>> import pyspeckit
        >>> x = pyspeckit.units.SpectroscopicAxis(np.linspace(-10,10,50), unit='km/s')
        >>> e = np.random.randn(50)
        >>> d = np.exp(-np.asarray(x)**2/2.)*5 + e
        >>> sp = pyspeckit.Spectrum(data=d, xarr=x, error=np.ones(50)*e.std())
        >>> sp.specfit(fittype='gaussian')
        >>> nwalkers = sp.specfit.fitter.npars * 2
        >>> emcee_ensemble = sp.specfit.fitter.get_emcee_ensemblesampler(sp.xarr, sp.data, sp.error, nwalkers)
        >>> p0 = np.array([sp.specfit.parinfo.values] * nwalkers)
        >>> p0 *= np.random.randn(*p0.shape) / 10. + 1.0
        >>> pos,logprob,state = emcee_ensemble.run_mcmc(p0,100)
        """
        try:
            import emcee
        except ImportError:
            return

        def probfunc(pars):
            return self.logp(xarr, data, error, pars=pars)

        sampler = emcee.EnsembleSampler(nwalkers,
                                        self.npars*self.npeaks+self.vheight,
                                        probfunc, **kwargs)

        return sampler

    def get_pymc(self, xarr, data, error, use_fitted_values=False, inf=np.inf,
                 use_adaptive=False, return_dict=False, **kwargs):
        """
        Create a pymc MCMC sampler.  Defaults to 'uninformative' priors

        Parameters
        ----------
        data : np.ndarray
        error : np.ndarray
        use_fitted_values : bool
            Each parameter with a measured error will have a prior defined by
            the Normal distribution with sigma = par.error and mean = par.value
        use_adaptive : bool
            Use the Adaptive Metropolis-Hastings sampler?

        Examples
        --------

        >>> x = pyspeckit.units.SpectroscopicAxis(np.linspace(-10,10,50), unit='km/s')
        >>> e = np.random.randn(50)
        >>> d = np.exp(-np.asarray(x)**2/2.)*5 + e
        >>> sp = pyspeckit.Spectrum(data=d, xarr=x, error=np.ones(50)*e.std())
        >>> sp.specfit(fittype='gaussian')
        >>> MCuninformed = sp.specfit.fitter.get_pymc(sp.xarr, sp.data, sp.error)
        >>> MCwithpriors = sp.specfit.fitter.get_pymc(sp.xarr, sp.data, sp.error, use_fitted_values=True)
        >>> MCuninformed.sample(1000)
        >>> MCuninformed.stats()['AMPLITUDE0']
        >>> # WARNING: This will fail because width cannot be set <0, but it may randomly reach that...
        >>> # How do you define a likelihood distribution with a lower limit?!
        >>> MCwithpriors.sample(1000)
        >>> MCwithpriors.stats()['AMPLITUDE0']

        """
        old_errsettings = np.geterr()
        try:
            import pymc
        finally:
            # pymc breaks error settings
            np.seterr(**old_errsettings)

        #def lowerlimit_like(x,lolim):
        #    "lower limit (log likelihood - set very positive for unacceptable values)"
        #    return (x>=lolim) / 1e10
        #def upperlimit_like(x,uplim):
        #    "upper limit"
        #    return (x<=uplim) / 1e10
        #LoLim = pymc.distributions.stochastic_from_dist('lolim', logp=lowerlimit_like, dtype=np.float, mv=False)
        #UpLim = pymc.distributions.stochastic_from_dist('uplim', logp=upperlimit_like, dtype=np.float, mv=False)

        funcdict = {}
        # very, very worrisome: pymc changes the values of parinfo
        parcopy = copy.deepcopy(self.parinfo)
        for par in parcopy:
            lolim = par.limits[0] if par.limited[0] else -inf
            uplim = par.limits[1] if par.limited[1] else inf
            if par.fixed:
                funcdict[par.parname] = pymc.distributions.Uniform(par.parname, par.value, par.value, value=par.value)
            elif use_fitted_values:
                if par.error > 0:
                    if any(par.limited):
                        try:
                            funcdict[par.parname] = pymc.distributions.TruncatedNormal(par.parname, par.value, 1./par.error**2, lolim, uplim)
                        except AttributeError:
                            # old versions used this?
                            funcdict[par.parname] = pymc.distributions.TruncNorm(par.parname, par.value, 1./par.error**2, lolim, uplim)
                    else:
                        funcdict[par.parname] = pymc.distributions.Normal(par.parname, par.value, 1./par.error**2)
                else:
                    if any(par.limited):
                        funcdict[par.parname] = pymc.distributions.Uniform(par.parname, lolim, uplim, value=par.value)
                    else:
                        funcdict[par.parname] = pymc.distributions.Uninformative(par.parname, value=par.value)
            elif any(par.limited):
                lolim = par.limits[0] if par.limited[0] else -1e10
                uplim = par.limits[1] if par.limited[1] else 1e10
                funcdict[par.parname] = pymc.distributions.Uniform(par.parname, lower=lolim, upper=uplim, value=par.value)
            else:
                funcdict[par.parname] = pymc.distributions.Uninformative(par.parname, value=par.value)

        d = dict(funcdict)

        def modelfunc(xarr, pars=parcopy, **kwargs):
            for k,v in kwargs.items():
                if k in list(pars.keys()):
                    pars[k].value = v

            return self.n_modelfunc(pars, **self.modelfunc_kwargs)(xarr)

        funcdict['xarr'] = xarr
        funcdet=pymc.Deterministic(name='f',eval=modelfunc,parents=funcdict,doc="The model function")
        d['f'] = funcdet

        datamodel = pymc.distributions.Normal('data', mu=funcdet,
                                              tau=1/np.asarray(error)**2,
                                              observed=True,
                                              value=np.asarray(data))
        d['data']=datamodel

        if return_dict:
            return d

        mc = pymc.MCMC(d)
        if use_adaptive:
            mc.use_step_method(pymc.AdaptiveMetropolis,[d[p] for p in self.parinfo.names])

        return mc

    def parse_3par_guesses(self, guesses):
        """
        Try to convert a set of interactive guesses (peak, center, width) into
        guesses appropriate to the model.
        """
        if len(guesses) % 3 != 0:
            raise ValueError("Guesses passed to parse_3par_guesses must have "
                             "length % 3 == 0")

        npeaks_guessed = len(guesses) // 3


        gtypes = [parse_offset_guess(gtype, gval)[0]
                  for gtype, gval in zip(itertools.cycle(self.guess_types),
                                         [0]*len(self.guess_types))]


        guess_dict = {(valid_guess_types[ii % 3], ii // 3): gval
                      for ii, gval in enumerate(guesses)}

        new_guesses = [guess_dict[(gtype, ii)]
                       if isinstance(gtype, str)
                       else gtype
                       for ii in range(npeaks_guessed)
                       for gtype in gtypes
                      ]

        new_guesses = [parse_offset_guess(gtype, gval)[1]
                       for gtype, gval in zip(itertools.cycle(self.guess_types),
                                              new_guesses)]

        assert len(new_guesses) % len(self.guess_types) == 0

        return new_guesses

class AstropyModel(SpectralModel):

    def __init__(self, model, shortvarnames=None, **kwargs):
        """
        Override the SpectralModel initialization
        """
        if hasattr(self,__doc__): # how do you extend a docstring really?
            self.__doc__ += SpectralModel.__doc__

        if shortvarnames is None:
            shortvarnames = model.param_names

        super(AstropyModel,self).__init__(model, len(model.parameters),
                                          shortvarnames=shortvarnames,
                                          model=model, **kwargs)

        self.mp = None
        self.vheight = False
        self.npeaks = 1


    def _make_parinfo(self, model=None):

        self.parinfo = ParinfoList([
            Parinfo(parname=name,value=value)
            for name,value in zip(model.param_names,model.parameters)])

        return self.parinfo, {}

    def _parse_parinfo(self, parinfo):
        """
        Parse a ParinfoList into astropy.models parameters
        """
        if len(parinfo) > self.npars:
            if len(parinfo) % self.npars != 0:
                raise ValueError("Need to have an integer number of models")
            else:
                self.modelfunc.param_names = parinfo.names
                self.modelfunc.parameters = parinfo.values
        else:
            self.modelfunc.param_names = parinfo.names
            self.modelfunc.parameters = parinfo.values

    def fitter(self, xax, data, err=None, quiet=True, veryverbose=False,
            debug=False, parinfo=None, params=None, npeaks=None, **kwargs):

        import astropy.models as models

        if npeaks is not None and npeaks > 1:
            raise NotImplementedError("Astropy models cannot be used to fit multiple peaks yet")

        if parinfo is not None:
            self._parse_parinfo(parinfo)
        if params is not None:
            self.modelfunc.parameters = params

        self.astropy_fitter = models.fitting.NonLinearLSQFitter(self.modelfunc)

        if err is None:
            self.astropy_fitter(xax, data, **kwargs)
        else:
            self.astropy_fitter(xax, data, weights=1./err**2, **kwargs)

        mpp = self.astropy_fitter.fitpars
        cov = self.astropy_fitter.covar
        if cov is None:
            mpperr = np.zeros(len(mpp))
        else:
            mpperr = cov.diagonal()
        self.model = self.astropy_fitter.model(xax)
        if err is None:
            chi2 = ((data-self.model)**2).sum()
        else:
            chi2 = ((data-self.model)**2/err**2).sum()

        # update object paramters
        self.modelfunc.parameters = mpp
        self._make_parinfo(self.modelfunc)

        return mpp,self.model,mpperr,chi2

    def n_modelfunc(self, pars=None, debug=False, **kwargs):
        """
        Only deals with single-peak functions
        """
        try:
            self._parse_parinfo(pars)
        except AttributeError:
            self.modelfunc.parameters = pars
        return self.modelfunc

def parse_offset_guess(gname, gval):
    """
    Utility function for handling guesses.  Allows guess types to be specified
    as 'amplitude*2' or 'width+3'.
    """
    operators = '+-*/'
    if not isinstance(gname, six.string_types):
        return gname, gval
    ops = [x for x in operators if x in gname]
    if len(ops)>1:
        raise ValueError("Invalid offset guess")
    elif len(ops) == 0:
        return gname,gval
    else:
        opmap = {"+": operator.add,
                 "-": operator.sub,
                 "*": operator.mul,
                 "/": operator.truediv,
                }
        op = ops[0]
        pars = gname.split(op)
        gname = [p for p in gname.split(op) if p in valid_guess_types][0]
        pars = [gval if p in valid_guess_types else float(p)
                for p in pars]
        gval = opmap[op](*pars)
        return gname, gval
