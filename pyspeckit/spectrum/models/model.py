"""
=============================
Generic SpectralModel wrapper 
=============================
.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>
"""
import numpy as np
from pyspeckit.mpfit import mpfit,mpfitException
import copy
import matplotlib.cbook as mpcb
import fitter
from . import mpfit_messages

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
            shortvarnames=("A","\\Delta x","\\sigma"), multisingle='multi', **kwargs):
        """
        modelfunc: the model function to be fitted.  Should take an X-axis (spectroscopic axis)
        as an input, followed by input parameters.
        npars - number of parameters required by the model

        parnames - a list or tuple of the parameter names

        parvalues - the initial guesses for the input parameters (defaults to ZEROS)

        parlimits - the upper/lower limits for each variable     (defaults to ZEROS)

        parfixed  - Can declare any variables to be fixed        (defaults to ZEROS)

        parerror  - technically an output parameter... hmm       (defaults to ZEROS)

        partied   - not the past tense of party.  Can declare, via text, that
            some parameters are tied to each other.  Defaults to zeros like the
            others, but it's not clear if that's a sensible default

        fitunits - convert X-axis to these units before passing to model

        parsteps - minimum step size for each paremeter          (defaults to ZEROS)

        npeaks   - default number of peaks to assume when fitting (can be overridden)

        shortvarnames - TeX names of the variables to use when annotating

        multisingle - Are there multiple peaks (no background will be fit) or
            just a single peak (a background may/will be fit)
        """

        self.modelfunc = modelfunc
        self.npars = npars 
        self.default_npars = npars
        self.multisingle = multisingle
        
        # this needs to be set once only
        self.shortvarnames = shortvarnames
        
        self.default_parinfo = None
        self.default_parinfo, kwargs = self._make_parinfo(**kwargs)
        self.parinfo = copy.copy(self.default_parinfo)

        self.modelfunc_kwargs = kwargs
        
    def _make_parinfo(self, params=None, parnames=None, parvalues=None,
            parlimits=None, parlimited=None, parfixed=None, parerror=None,
            partied=None, fitunits=None, parsteps=None, npeaks=1,
            parinfo=None,
            names=None, values=None, limits=None,
            limited=None, fixed=None, error=None, tied=None, steps=None,
            negamp=None,
            limitedmin=None, limitedmax=None,
            minpars=None, maxpars=None,
            vheight=False,
            debug=False,
            **kwargs):

        # for backwards compatibility - partied = tied, etc.
        for varname in str.split("parnames,parvalues,parsteps,parlimits,parlimited,parfixed,parerror,partied",","):
            shortvarname = varname.replace("par","")
            if locals()[shortvarname] is not None:
                # HACK!  locals() failed for unclear reasons...
                exec("%s = %s" % (varname,shortvarname))

        if params is not None and parvalues is not None:
            raise ValueError("parvalues and params both specified; they're redundant so that's not allowed.")
        elif params is not None and parvalues is None:
            parvalues = params

        if parnames is not None: 
            self.parnames = parnames
        elif parnames is None and self.parnames is not None:
            parnames = self.parnames
        elif self.default_parinfo is not None:
            parnames = [p['parname'] for p in self.default_parinfo]

        if limitedmin is not None:
            if limitedmax is not None:
                parlimited = zip(limitedmin,limitedmax)
            else:
                parlimited = zip(limitedmin,(False,)*len(parnames))
        elif limitedmax is not None:
            parlimited = zip((False,)*len(parnames),limitedmax)
        elif self.default_parinfo is not None:
            parlimited = [p['limited'] for p in self.default_parinfo]

        if minpars is not None:
            if maxpars is not None:
                parlimits = zip(minpars,maxpars)
            else:
                parlimits = zip(minpars,(False,)*len(parnames))
        elif maxpars is not None:
            parlimits = zip((False,)*len(parnames),maxpars)
        elif self.default_parinfo is not None:
            parlimits = [p['limits'] for p in self.default_parinfo]

        self.fitunits = fitunits
        self.npeaks = npeaks

        # the height / parvalue popping needs to be done before the temp_pardict is set in order to make sure
        # that the height guess isn't assigned to the amplitude
        self.vheight = vheight
        if vheight and len(self.parinfo) == self.default_npars and len(parvalues) == self.default_npars + 1:
            # if the right number of parameters are passed, the first is the height
            self.parinfo = [ {'n':0, 'value':parvalues.pop(0), 'limits':(0,0),
                'limited': (False,False), 'fixed':False, 'parname':'HEIGHT',
                'error': 0, 'tied':"" } ]
        elif vheight and len(self.parinfo) == self.default_npars and len(parvalues) == self.default_npars:
            # if you're one par short, guess zero
            self.parinfo = [ {'n':0, 'value': 0, 'limits':(0,0),
                'limited': (False,False), 'fixed':False, 'parname':'HEIGHT',
                'error': 0, 'tied':"" } ]
        elif vheight and len(self.parinfo) == self.default_npars+1 and len(parvalues) == self.default_npars+1:
            # the right numbers are passed *AND* there is already a height param
            self.parinfo = [ {'n':0, 'value':parvalues.pop(0), 'limits':(0,0),
                'limited': (False,False), 'fixed':False, 'parname':'HEIGHT',
                'error': 0, 'tied':"" } ]
            #heightparnum = (i for i,s in self.parinfo if 'HEIGHT' in s['parname'])
            #for hpn in heightparnum:
            #    self.parinfo[hpn]['value'] = parvalues[0]
        elif vheight:
            raise ValueError('VHEIGHT is specified but a case was found that did not allow it to be included.')
        else:
            self.parinfo = []

        if debug: print "After VHEIGHT parse len(parinfo): %i   vheight: %s" % (len(self.parinfo), vheight)


        # this is a clever way to turn the parameter lists into a dict of lists
        # clever = hard to read
        temp_pardict = dict([(varname, np.zeros(self.npars*self.npeaks, dtype='bool'))
            if locals()[varname] is None else (varname, list(locals()[varname]) )
            for varname in str.split("parnames,parvalues,parsteps,parlimits,parlimited,parfixed,parerror,partied",",")])
        temp_pardict['parlimits'] = parlimits if parlimits is not None else [(0,0)] * (self.npars*self.npeaks)
        temp_pardict['parlimited'] = parlimited if parlimited is not None else [(False,False)] * (self.npars*self.npeaks)
        for k,v in temp_pardict.iteritems():
            if (self.npars*self.npeaks) / len(v) > 1:
                temp_pardict[k] = list(v) * ((self.npars*self.npeaks) / len(v))

        # generate the parinfo dict
        # note that 'tied' must be a blank string (i.e. ""), not False, if it is not set
        # parlimited, parfixed, and parlimits are all two-element items (tuples or lists)
        self.parinfo += [ {'n':ii+self.npars*jj+vheight,
            'value':temp_pardict['parvalues'][ii+self.npars*jj],
            'step':temp_pardict['parsteps'][ii+self.npars*jj],
            'limits':temp_pardict['parlimits'][ii+self.npars*jj],
            'limited':temp_pardict['parlimited'][ii+self.npars*jj],
            'fixed':temp_pardict['parfixed'][ii+self.npars*jj],
            'parname':temp_pardict['parnames'][ii].upper()+"%0i" % jj,
            'error':temp_pardict['parerror'][ii+self.npars*jj],
            'tied':temp_pardict['partied'][ii+self.npars*jj] if temp_pardict['partied'][ii+self.npars*jj] else ""} 
            for jj in xrange(self.npeaks)
            for ii in xrange(self.npars) ] # order matters!

        if debug: print "After Generation step len(parinfo): %i   vheight: %s" % (len(self.parinfo), vheight)

        if debug > True: import pdb; pdb.set_trace()

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

        return self.parinfo, kwargs

    def n_modelfunc(self, pars, **kwargs):
        """
        Simple wrapper to deal with N independent peaks for a given spectral model
        """
        pars = list(pars)
        def L(x):
            v = np.zeros(len(x))
            if self.vheight: v += pars[0]
            # use len(pars) instead of self.npeaks because we want this to work
            # independent of the current best fit
            for jj in xrange((len(pars)-self.vheight)/self.npars):
                v += self.modelfunc(x, *pars[jj*self.npars+self.vheight:(jj+1)*self.npars+self.vheight], **kwargs)
            return v
        return L

    def mpfitfun(self,x,y,err=None):
        """
        Wrapper function to compute the fit residuals in an mpfit-friendly format
        """
        if err is None:
            def f(p,fjac=None): return [0,(y-self.n_modelfunc(p, **self.modelfunc_kwargs)(x))]
        else:
            def f(p,fjac=None): return [0,(y-self.n_modelfunc(p, **self.modelfunc_kwargs)(x))/err]
        return f

    def __call__(self, *args, **kwargs):
        if self.multisingle == 'single':
            # Generate a variable-height version of the model
            func = fitter.vheightmodel(self.modelfunc)
            return self.fitter(*args, **kwargs)
        elif self.multisingle == 'multi':
            return self.fitter(*args,**kwargs)


    def fitter(self, xax, data, err=None, quiet=True, veryverbose=False,
            debug=False, parinfo=None, **kwargs):
        """
        Run the fitter.  Must pass the x-axis and data.  Can include
        error, parameter guesses, and a number of verbosity parameters.

        quiet - pass to mpfit.  If False, will print out the parameter values
            for each iteration of the fitter

        veryverbose - print out a variety of mpfit output parameters

        debug - raise an exception (rather than a warning) if chi^2 is nan

        accepts *tied*, *limits*, *limited*, and *fixed* as keyword arguments.
            They must be lists of length len(params)

        parinfo - You can override the class parinfo dict with this, though
            that largely defeats the point of having the wrapper class.  This class
            does NO checking for whether the parinfo dict is valid.

        kwargs are passed to mpfit after going through _make_parinfo to strip out things
        used by this class
        """

        if parinfo is None:
            parinfo, kwargs = self._make_parinfo(debug=debug, **kwargs)
        else:
            if debug: print "Using user-specified parinfo dict"
            # clean out disallowed kwargs (don't want to pass them to mpfit)
            throwaway, kwargs = self._make_parinfo(debug=debug, **kwargs)

        self.xax = xax # the 'stored' xax is just a link to the original
        if hasattr(xax,'convert_to_unit') and self.fitunits is not None:
            # some models will depend on the input units.  For these, pass in an X-axis in those units
            # (gaussian, voigt, lorentz profiles should not depend on units.  Ammonia, formaldehyde,
            # H-alpha, etc. should)
            xax = copy.copy(xax)
            xax.convert_to_unit(self.fitunits, quiet=quiet)

        if np.any(np.isnan(data)) or np.any(np.isinf(data)):
            err[np.isnan(data) + np.isinf(data)] = np.inf
            data[np.isnan(data) + np.isinf(data)] = 0

        if debug:
            for p in parinfo: print p
            print "\n".join(["%s %i: tied: %s value: %s" % (p['parname'],p['n'],p['tied'],p['value']) for p in parinfo])

        mp = mpfit(self.mpfitfun(xax,data,err),parinfo=parinfo,quiet=quiet,**kwargs)
        mpp = mp.params
        if mp.perror is not None: mpperr = mp.perror
        else: mpperr = mpp*0
        chi2 = mp.fnorm

        if mp.status == 0:
            raise mpfitException(mp.errmsg)

        for i,(p,e) in enumerate(zip(mpp,mpperr)):
            self.parinfo[i]['value'] = p
            self.parinfo[i]['error'] = e

        if veryverbose:
            print "Fit status: ",mp.status
            print "Fit error message: ",mp.errmsg
            print "Fit message: ",mpfit_messages[mp.status]
            for i,p in enumerate(mpp):
                print self.parinfo[i]['parname'],p," +/- ",mpperr[i]
            print "Chi2: ",mp.fnorm," Reduced Chi2: ",mp.fnorm/len(data)," DOF:",len(data)-len(mpp)

        self.mp = mp
        self.mpp = mpp
        self.mpperr = mpperr
        self.model = self.n_modelfunc(mpp,**self.modelfunc_kwargs)(xax)
        if debug:
            print "Modelpars: ",self.mpp
        if np.isnan(chi2):
            if debug:
                raise ValueError("Error: chi^2 is nan")
            else:
                print "Warning: chi^2 is nan"
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

    def annotations(self, shortvarnames=None):
        """
        Return a list of TeX-formatted labels
        """
        from decimal import Decimal # for formatting
        svn = self.shortvarnames if shortvarnames is None else shortvarnames
        # if pars need to be replicated....
        if len(svn) < self.npeaks*self.npars:
            svn = svn * self.npeaks
        label_list = [(
                "$%s(%i)$=%8s $\\pm$ %8s" % (svn[ii+jj*self.npars],jj,
                Decimal("%g" % self.mpp[ii+jj*self.npars+self.vheight]).quantize(Decimal("%0.2g" % (min(self.mpp[ii+jj*self.npars+self.vheight],self.mpperr[ii+jj*self.npars+self.vheight])))),
                Decimal("%g" % self.mpperr[ii+jj*self.npars+self.vheight]).quantize(Decimal("%0.2g" % (self.mpperr[ii+jj*self.npars+self.vheight]))),)
                          ) for jj in range(self.npeaks) for ii in range(self.npars)]
        labels = tuple(mpcb.flatten(label_list))
        return labels

    def components(self, xarr, pars, **kwargs):
        """
        Return a numpy ndarray of the independent components of the fits
        """

        modelcomponents = np.array(
            [self.modelfunc(xarr,
                *pars[i*self.npars:(i+1)*self.npars],
                return_components=True,
                **self.modelfunc_kwargs)
            for i in range(self.npeaks)])

        return modelcomponents

    def integral(self, modelpars, **kwargs):
        """
        Extremely simple integrator:
        IGNORES modelpars;
        just sums self.model
        """

        return self.model.sum()
