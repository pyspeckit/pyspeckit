"""
=============================
Generic SpectralModel wrapper 
=============================
.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>
"""
import numpy as np
from mpfit import mpfit
import copy
import matplotlib.cbook as mpcb
import fitter

class SpectralModel(fitter.SimpleFitter):
    """
    A wrapper class for a spectra model.  Includes internal functions to
    generate multi-component models, annotations, integrals, and individual
    components.  The declaration can be complex, since you should name
    individual variables, set limits on them, set the units the fit will be
    performed in, and set the annotations to be used.  Check out some
    of the hyperfine codes (hcn, n2hp) for examples.
    """

    def __init__(self, modelfunc, npars, parnames=None, parvalues=None,
            parlimits=None, parlimited=None, parfixed=None, parerror=None,
            partied=None, fitunits=None, parsteps=None, npeaks=1,
            shortvarnames=("A","v","\\sigma"), parinfo=None, **kwargs):
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
        """

        self.modelfunc = modelfunc
        self.npars = npars
        self.parnames = parnames
        self.fitunits = fitunits
        self.npeaks = npeaks
        self.shortvarnames = shortvarnames

        temp_pardict = dict([(varname, np.zeros(self.npars, dtype='bool')) if locals()[varname] is None else (varname, locals()[varname])
            for varname in str.split("parnames,parvalues,parsteps,parlimits,parlimited,parfixed,parerror,partied",",")])

        if parinfo is not None:
            self.parinfo = parinfo
        else:
            # generate the parinfo dict
            # note that 'tied' must be a blank string (i.e. ""), not False, if it is not set
            # parlimited, parfixed, and parlimits are all two-element items (tuples or lists)
            self.parinfo = [ {'n':ii,
                'value':temp_pardict['parvalues'][ii],
                'step':temp_pardict['parsteps'][ii],
                'limits':temp_pardict['parlimits'][ii],
                'limited':temp_pardict['parlimited'][ii],
                'fixed':temp_pardict['parfixed'][ii],
                'parname':temp_pardict['parnames'][ii]+"%0i" % jj,
                'error':temp_pardict['parerror'][ii],
                'tied':temp_pardict['partied'][ii] if temp_pardict['partied'][ii] else ""} 
                for ii in xrange(self.npars) for jj in xrange(self.npeaks)]

        self.modelfunc_kwargs = kwargs

    def n_modelfunc(self, pars, **kwargs):
        """
        Simple wrapper to deal with N independent peaks for a given spectral model
        """
        def L(x):
            v = np.zeros(len(x))
            for jj in xrange(self.npeaks):
                v += self.modelfunc(x, *pars[jj*self.npars:(jj+1)*self.npars], **kwargs)
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

    def __call__(self, xax, data, err=None, params=(), quiet=True, 
            veryverbose=False, npeaks=None, debug=False, **kwargs):
        """
        Run the fitter.  Must pass the x-axis and data.  Can include
        error, parameter guesses, and a number of verbosity parameters.

        quiet - pass to mpfit.  If False, will print out the parameter values
            for each iteration of the fitter

        veryverbose - print out a variety of mpfit output parameters

        debug - raise an exception (rather than a warning) if chi^2 is nan

        kwargs are passed to mpfit
        """

        if npeaks is not None:
            if npeaks > self.npeaks:
                # duplicate the current parameters npeaks-oldnpeaks times
                for ii in xrange(npeaks-self.npeaks):
                    self.parinfo += self.parinfo[:self.npars]
                    self.shortvarnames += self.shortvarnames
                    # replace the number for each parameter
                    for jj in xrange(self.npars): 
                        self.parinfo[self.npars*(ii+1)+jj]['n'] = self.parinfo[jj]['n'] + self.npars
                        self.parinfo[self.npars*(ii+1)+jj]['parname'] = self.parinfo[jj]['parname'].replace('0','%0i' % (ii+1))
            self.npeaks = npeaks

        if len(params) == self.npars*self.npeaks:
            for par,guess in zip(self.parinfo,params):
                par['value'] = guess
        
        for varname in str.split("limits,limited,fixed,tied",","):
            if varname in kwargs:
                var = kwargs.pop(varname)
                for pi in self.parinfo:
                    pi[varname] = var[pi['n']]

        if hasattr(xax,'convert_to_unit') and self.fitunits is not None:
            # some models will depend on the input units.  For these, pass in an X-axis in those units
            # (gaussian, voigt, lorentz profiles should not depend on units.  Ammonia, formaldehyde,
            # H-alpha, etc. should)
            xax = copy.copy(xax)
            xax.convert_to_unit(self.fitunits, quiet=quiet)

        if np.any(np.isnan(data)) or np.any(np.isinf(data)):
            err[np.isnan(data) + np.isinf(data)] = np.inf
            data[np.isnan(data) + np.isinf(data)] = 0

        mp = mpfit(self.mpfitfun(xax,data,err),parinfo=self.parinfo,quiet=quiet,**kwargs)
        mpp = mp.params
        if mp.perror is not None: mpperr = mp.perror
        else: mpperr = mpp*0
        chi2 = mp.fnorm

        if mp.status == 0:
            raise Exception(mp.errmsg)

        if veryverbose:
            print "Fit status: ",mp.status
            print "Fit error message: ",mp.errmsg
            print "Fit message: ",mpfit_messages[mp.status]
            for i,p in enumerate(mpp):
                self.parinfo[i]['value'] = p
                print self.parinfo[i]['parname'],p," +/- ",mpperr[i]
            print "Chi2: ",mp.fnorm," Reduced Chi2: ",mp.fnorm/len(data)," DOF:",len(data)-len(mpp)

        self.mp = mp
        self.mpp = mpp#[1:]
        self.mpperr = mpperr#[1:]
        self.model = self.n_modelfunc(mpp,**self.modelfunc_kwargs)(xax)
        if np.isnan(chi2):
            if debug:
                raise ValueError("Error: chi^2 is nan")
            else:
                print "Warning: chi^2 is nan"
        return mpp,self.model,mpperr,chi2


    def annotations(self, shortvarnames=None):
        """
        Return a list of TeX-formatted labels
        """
        from decimal import Decimal # for formatting
        svn = self.shortvarnames if shortvarnames is None else shortvarnames
        label_list = [(
                "$%s(%i)$=%8s $\\pm$ %8s" % (svn[ii+jj*self.npars],jj,
                Decimal("%g" % self.mpp[ii+jj*self.npars]).quantize(Decimal("%0.2g" % (min(self.mpp[ii+jj*self.npars],self.mpperr[ii+jj*self.npars])))),
                Decimal("%g" % self.mpperr[ii+jj*self.npars]).quantize(Decimal("%0.2g" % (self.mpperr[ii+jj*self.npars]))),)
                          ) for jj in range(self.npeaks) for ii in range(self.npars)]
        labels = tuple(mpcb.flatten(label_list))
        return labels

    def components(self, xarr, pars):
        """
        Return a numpy ndarray of the independent components of the fits
        """

        modelcomponents = np.concatenate(
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
