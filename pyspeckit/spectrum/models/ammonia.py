"""
========================================
Ammonia inversion transition TROT fitter
========================================

Ammonia inversion transition TROT fitter translated from Erik Rosolowsky's
https://github.com/low-sky/nh3fit

.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>

Module API
^^^^^^^^^^

"""
from __future__ import division

import numpy as np
from ...mpfit import mpfit
from ...spectrum.parinfo import ParinfoList,Parinfo
from . import fitter
from . import model
import matplotlib.cbook as mpcb
import copy
from astropy import log
from six import iteritems
from . import mpfit_messages
import operator
import string
import warnings

from .ammonia_constants import (line_names, freq_dict, aval_dict, ortho_dict,
                                TCMB,
                                voff_lines_dict, tau_wts_dict)


def ammonia(xarr, trot=20, tex=None, ntot=14, width=1, xoff_v=0.0, fortho=0.0,
            tau=None, fillingfraction=None, return_tau=False,
            return_tau_profile=False, background_tb=TCMB, verbose=False,
            return_components=False, debug=False, line_names=line_names,
            ignore_neg_models=False):
    """
    Generate a model Ammonia spectrum based on input temperatures, column, and
    gaussian parameters.  The returned model will be in Kelvin (brightness
    temperature) units.

    Note that astropy units are not used internally for performance reasons.  A
    wrapped version of this module including those units would be a good idea,
    as it is definitely possible to implement this with unit support and good
    performance.

    Parameters
    ----------
    xarr: `pyspeckit.spectrum.units.SpectroscopicAxis`
        Array of wavelength/frequency values
    trot: float
        The rotational temperature of the lines.  This is the excitation
        temperature that governs the relative populations of the rotational
        states.
    tex: float or None
        Excitation temperature. Assumed LTE if unspecified (``None``) or if
        tex>trot.  This is the excitation temperature for *all* of the modeled
        lines, which means we are explicitly assuming T_ex is the same for all
        lines.
    ntot: float
        Total log column density of NH3.  Can be specified as a float in the
        range 5-25
    width: float
        Line width (Gaussian sigma) in km/s
    xoff_v: float
        Line offset in km/s
    fortho: float
        Fraction of NH3 molecules in ortho state.  Default assumes all para
        (fortho=0).
    tau: None or float
        If tau (optical depth in the 1-1 line) is specified, ntot is NOT fit
        but is set to a fixed value.  The optical depths of the other lines are
        fixed relative to tau_oneone
    fillingfraction: None or float
        fillingfraction is an arbitrary scaling factor to apply to the model
    return_tau: bool
        Return a dictionary of the optical depths in each line instead of a
        synthetic spectrum
    return_tau_profile: bool
        Return a dictionary of the optical depth profiles in each line, i.e.,
        the optical depths that will be used in conjunction with T_ex to produce
        the synthetic spectrum
    return_components: bool
        Return a list of arrays, one for each hyperfine component, instead of
        just one array
    background_tb : float
        The background brightness temperature.  Defaults to TCMB.
    ignore_neg_models: bool
        Normally if background=TCMB and the model is negative, an exception
        will be raised.  This parameter will simply skip that exception.  Use
        with extreme caution: negative models (absorption spectra against the
        CMB) are not physical!  You may want to allow this in some cases
        because there can be numerical issues where the model goes negative
        when it shouldn't.
    verbose: bool
        More messages
    debug: bool
        For debugging.

    Returns
    -------
    spectrum: `numpy.ndarray`
        Synthetic spectrum with same shape as ``xarr``
    component_list: list
        List of `numpy.ndarray`'s, one for each hyperfine component
    tau_dict: dict
        Dictionary of optical depth values for the various lines
        (if ``return_tau`` is set)
    """

    from .ammonia_constants import (ckms, ccms, h, kb,
                                    Jortho, Jpara, Brot, Crot)

    # Convert X-units to frequency in GHz
    if xarr.unit.to_string() != 'GHz':
        xarr = xarr.as_unit('GHz')

    if tex is None:
        log.warning("Assuming tex=trot")
        tex = trot
    elif isinstance(tex, dict):
        for k in tex:
            assert k in line_names,"{0} not in line list".format(k)
        line_names = tex.keys()
    elif tex > trot:
        warnings.warn("tex > trot in the ammonia model.  "
                      "This is unphysical and "
                      "suggests that you may need to constrain tex.  See "
                      "ammonia_model_restricted_tex.")
    if width < 0:
        return np.zeros(xarr.size)*np.nan
    elif width == 0:
        return np.zeros(xarr.size)

    from .ammonia_constants import line_name_indices, line_names as original_line_names

    # recreate line_names keeping only lines with a specified tex
    # using this loop instead of tex.keys() preserves the order & data type
    line_names = [k for k in original_line_names if k in line_names]

    if 5 <= ntot <= 25:
        # allow ntot to be specified as a logarithm.  This is
        # safe because ntot < 1e10 gives a spectrum of all zeros, and the
        # plausible range of columns is not outside the specified range
        lin_ntot = 10**ntot
    else:
        raise ValueError("ntot, the logarithmic total column density,"
                         " must be in the range 5 - 25")

    tau_dict = {}

    """
    Column density is the free parameter.  It is used in conjunction with
    the full partition function to compute the optical depth in each band
    """
    Zpara = (2*Jpara+1)*np.exp(-h*(Brot*Jpara*(Jpara+1)+
                                   (Crot-Brot)*Jpara**2)/(kb*trot))
    Zortho = 2*(2*Jortho+1)*np.exp(-h*(Brot*Jortho*(Jortho+1)+
                                       (Crot-Brot)*Jortho**2)/(kb*trot))
    Qpara = Zpara.sum()
    Qortho = Zortho.sum()

    log.debug("Partition Function: Q_ortho={0}, Q_para={1}".format(Qortho, Qpara))

    for linename in line_names:
        if ortho_dict[linename]:
            # define variable "ortho_or_para_frac" that will be the ortho
            # fraction in the case of an ortho transition or the para
            # fraction for a para transition
            ortho_or_parafrac = fortho
            Z = Zortho
            Qtot = Qortho
        else:
            ortho_or_parafrac = 1.0-fortho
            Z = Zpara
            Qtot = Qpara

        # for a complete discussion of these equations, please see
        # https://github.com/keflavich/pyspeckit/blob/ammonia_equations/examples/AmmoniaLevelPopulation.ipynb
        # https://github.com/pyspeckit/pyspeckit/blob/master/examples/AmmoniaLevelPopulation.ipynb
        # and
        # http://low-sky.github.io/ammoniacolumn/
        # and
        # https://github.com/pyspeckit/pyspeckit/pull/136

        # short variable names for readability
        frq = freq_dict[linename]
        partition = Z[line_name_indices[linename]]
        aval = aval_dict[linename]

        # Total population of the higher energy inversion transition
        population_rotstate = lin_ntot * ortho_or_parafrac * partition/Qtot

        if isinstance(tex, dict):
            expterm = ((1-np.exp(-h*frq/(kb*tex[linename]))) /
                       (1+np.exp(-h*frq/(kb*tex[linename]))))
        else:
            expterm = ((1-np.exp(-h*frq/(kb*tex))) /
                       (1+np.exp(-h*frq/(kb*tex))))
        fracterm = (ccms**2 * aval / (8*np.pi*frq**2))
        widthterm = (ckms/(width*frq*(2*np.pi)**0.5))

        tau_i = population_rotstate * fracterm * expterm * widthterm
        tau_dict[linename] = tau_i

        log.debug("Line {0}: tau={1}, expterm={2}, pop={3},"
                  " partition={4}"
                  .format(linename, tau_i, expterm, population_rotstate,
                          partition))

    # allow tau(11) to be specified instead of ntot
    # in the thin case, this is not needed: ntot plays no role
    # this process allows you to specify tau without using the approximate equations specified
    # above.  It should remove ntot from the calculations anyway...
    if tau is not None:
        tau11_temp = tau_dict['oneone']
        # re-scale all optical depths so that tau is as specified, but the relative taus
        # are sest by the kinetic temperature and partition functions
        for linename,t in iteritems(tau_dict):
            tau_dict[linename] = t * tau/tau11_temp

    if return_tau:
        return tau_dict


    model_spectrum = _ammonia_spectrum(xarr, tex, tau_dict, width, xoff_v,
                                       fortho, line_names,
                                       background_tb=background_tb,
                                       fillingfraction=fillingfraction,
                                       return_components=return_components,
                                       return_tau_profile=return_tau_profile
                                      )

    if not return_tau_profile and model_spectrum.min() < 0 and background_tb == TCMB and not ignore_neg_models:
        raise ValueError("Model dropped below zero.  That is not possible "
                         " normally.  Here are the input values: "+
                         ("tex: {0} ".format(tex)) +
                         ("trot: %f " % trot) +
                         ("ntot: %f " % ntot) +
                         ("width: %f " % width) +
                         ("xoff_v: %f " % xoff_v) +
                         ("fortho: %f " % fortho)
                         )


    if verbose or debug:
        log.info("trot: %g  tex: %s  ntot: %g  width: %g  xoff_v: %g  "
                 "fortho: %g  fillingfraction: %g" % (trot, tex, ntot, width,
                                                      xoff_v, fortho,
                                                      fillingfraction))


    return model_spectrum

def cold_ammonia(xarr, tkin, **kwargs):
    """
    Generate a model Ammonia spectrum based on input temperatures, column, and
    gaussian parameters

    Parameters
    ----------
    xarr: `pyspeckit.spectrum.units.SpectroscopicAxis`
        Array of wavelength/frequency values
    tkin: float
        The kinetic temperature of the lines in K.  Will be converted to
        rotational temperature following the scheme of Swift et al 2005
        (http://esoads.eso.org/abs/2005ApJ...620..823S, eqn A6) and further
        discussed in Equation 7 of Rosolowsky et al 2008
        (http://adsabs.harvard.edu/abs/2008ApJS..175..509R)
    """

    dT0 = 41.18 # Energy difference between (2,2) and (1,1) in K
    trot = tkin * (1 + (tkin/dT0)*np.log(1 + 0.6*np.exp(-15.7/tkin)))**-1
    log.debug("Cold ammonia turned T_K = {0} into T_rot = {1}".format(tkin,trot))

    return ammonia(xarr, trot=trot, **kwargs)

def ammonia_thin(xarr, tkin=20, tex=None, ntot=14, width=1, xoff_v=0.0,
                 fortho=0.0, tau=None, return_tau=False, **kwargs):
    """
    Use optical depth in the 1-1 line as a free parameter
    The optical depths of the other lines are then set by the kinetic
    temperature

    tkin is used to compute trot assuming a 3-level system consisting of (1,1),
    (2,1), and (2,2) as in Swift et al, 2005 [2005ApJ...620..823S]
    """

    tau_dict = {}

    tex = tkin

    dT0 = 41.5                    # Energy diff between (2,2) and (1,1) in K
    trot = tkin/(1+tkin/dT0*np.log(1+0.6*np.exp(-15.7/tkin)))
    tau_dict['oneone'] = tau
    tau_dict['twotwo'] = tau*(23.722/23.694)**2*4/3.*5/3.*np.exp(-41.5/trot)
    tau_dict['threethree'] = tau*(23.8701279/23.694)**2*3/2.*14./3.*np.exp(-101.1/trot)
    tau_dict['fourfour'] = tau*(24.1394169/23.694)**2*8/5.*9/3.*np.exp(-177.34/trot)
    line_names = tau_dict.keys()
    # TODO: Raise a warning if tkin > (some value), probably 50 K, because
    # the 3-level system approximation used here will break down.

    if return_tau:
        return tau_dict
    else:
        return _ammonia_spectrum(xarr, tex, tau_dict, width, xoff_v, fortho,
                                 line_names, **kwargs)

def _ammonia_spectrum(xarr, tex, tau_dict, width, xoff_v, fortho, line_names,
                      background_tb=TCMB, fillingfraction=None,
                      return_components=False, return_tau_profile=False):
    """
    Helper function: given a dictionary of ammonia optical depths,
    an excitation tmeperature etc, produce the spectrum.

    The default return units are brightness temperature in Kelvin.  If
    ``return_tau_profile`` is specified, the returned "spectrum" will be
    a spectrum of optical depths, not an intensity spectrum.

    If ``return_components`` is specified, a list of spectra will be returned,
    where each spectrum represents one of the hyperfine components of the
    particular ammonia line being modeled.
    """
    from .ammonia_constants import (ckms, h, kb)

    # fillingfraction is an arbitrary scaling for the data
    # The model will be (normal model) * fillingfraction
    if fillingfraction is None:
        fillingfraction = 1.0

    # "runspec" means "running spectrum": it is accumulated over a loop
    runspec = np.zeros(len(xarr))

    if return_components:
        components = []

    if return_tau_profile:
        tau_profile = {}

    for linename in line_names:
        voff_lines = np.array(voff_lines_dict[linename])
        tau_wts = np.array(tau_wts_dict[linename])

        lines = (1-voff_lines/ckms)*freq_dict[linename]/1e9
        tau_wts = tau_wts / (tau_wts).sum()
        nuwidth = np.abs(width/ckms*lines)
        nuoff = xoff_v/ckms*lines

        # tau array
        tauprof = np.zeros(len(xarr))
        for kk,nuo in enumerate(nuoff):
            tauprof_ = (tau_dict[linename] * tau_wts[kk] *
                        np.exp(-(xarr.value+nuo-lines[kk])**2 /
                               (2.0*nuwidth[kk]**2)))
            if return_components:
                components.append(tauprof_)
            tauprof += tauprof_

        if return_tau_profile:
            tau_profile[linename] = tauprof

        T0 = (h*xarr.value*1e9/kb) # "temperature" of wavelength

        if isinstance(tex, dict):
            runspec = ((T0/(np.exp(T0/tex[linename])-1) -
                        T0/(np.exp(T0/background_tb)-1)) *
                       (1-np.exp(-tauprof)) * fillingfraction + runspec)
        else:
            runspec = ((T0/(np.exp(T0/tex)-1) -
                        T0/(np.exp(T0/background_tb)-1)) *
                       (1-np.exp(-tauprof)) * fillingfraction + runspec)


    if return_components:
        if isinstance(tex, dict):
            term1 = [(T0/(np.exp(T0/tex[linename])-1)-T0/(np.exp(T0/background_tb)-1))
                     for linename in line_names]
        else:
            term1 = (T0/(np.exp(T0/tex)-1)-T0/(np.exp(T0/background_tb)-1))
        return term1*(1-np.exp(-1*np.array(components)))
    elif return_tau_profile:
        return tau_profile
    else:
        return runspec



class ammonia_model(model.SpectralModel):
    """
    The basic Ammonia (NH3) model with 6 free parameters:
        Trot, Tex, ntot, width, xoff_v, and fortho

    Trot is the rotational temperature.  It governs the relative populations of
    the rotational states, i.e., the relative strength of different transitions

    Tex is the excitation temperature.  It is assumed constant across all
    states, which is not always a good assumption - a radiative transfer and
    excitation model is required to constrain this, though.

    ntot is the total column density of p-NH3 integrated over all states.

    width is the linewidth (the Gaussian sigma)

    xoff_v is the velocity offset / line of sight velocity

    fortho is the ortho fraction  (northo / (northo+npara))
    """

    def __init__(self,npeaks=1,npars=6,
                 parnames=['trot','tex','ntot','width','xoff_v','fortho'],
                 **kwargs):
        npeaks = self.npeaks = int(npeaks)
        npars = self.npars = int(npars)
        self._default_parnames = parnames
        self.parnames = copy.copy(self._default_parnames)

        # all fitters must have declared modelfuncs, which should take the fitted pars...
        self.modelfunc = ammonia
        self.n_modelfunc = self.n_ammonia

        # for fitting ammonia simultaneously with a flat background
        self.onepeakammonia = fitter.vheightmodel(ammonia)
        #self.onepeakammoniafit = self._fourparfitter(self.onepeakammonia)

        self.default_parinfo = None
        self.default_parinfo, kwargs = self._make_parinfo(**kwargs)

        # Remove keywords parsed by parinfo and ignored by the fitter
        for kw in ('tied','partied'):
            if kw in kwargs:
                kwargs.pop(kw)

        # enforce ammonia-specific parameter limits
        for par in self.default_parinfo:
            if 'tex' in par.parname.lower():
                par.limited = (True,par.limited[1])
                par.limits = (max(par.limits[0],TCMB), par.limits[1])
            if 'trot' in par.parname.lower():
                par.limited = (True,par.limited[1])
                par.limits = (max(par.limits[0],TCMB), par.limits[1])
            if 'width' in par.parname.lower():
                par.limited = (True,par.limited[1])
                par.limits = (max(par.limits[0],0), par.limits[1])
            if 'fortho' in par.parname.lower():
                par.limited = (True,True)
                if par.limits[1] != 0:
                    par.limits = (max(par.limits[0],0), min(par.limits[1],1))
                else:
                    par.limits = (max(par.limits[0],0), 1)
            if 'ntot' in par.parname.lower():
                par.limited = (True,par.limited[1])
                par.limits = (max(par.limits[0],0), par.limits[1])

        self.parinfo = copy.copy(self.default_parinfo)

        self.modelfunc_kwargs = kwargs
        # lower case? self.modelfunc_kwargs.update({'parnames':self.parinfo.parnames})
        self.use_lmfit = kwargs.pop('use_lmfit') if 'use_lmfit' in kwargs else False
        self.fitunit = 'GHz'

    def __call__(self,*args,**kwargs):
        return self.multinh3fit(*args,**kwargs)

    def n_ammonia(self, pars=None, parnames=None, **kwargs):
        """
        Returns a function that sums over N ammonia line profiles, where N is the length of
        trot,tex,ntot,width,xoff_v,fortho *OR* N = len(pars) / 6

        The background "height" is assumed to be zero (you must "baseline" your
        spectrum before fitting)

        *pars* [ list ]
            a list with len(pars) = (6-nfixed)n, assuming
            trot,tex,ntot,width,xoff_v,fortho repeated

        *parnames* [ list ]
            len(parnames) must = len(pars).  parnames determine how the ammonia
            function parses the arguments
        """
        npeaks = self.npeaks
        npars = len(self.default_parinfo)

        if hasattr(pars,'values'):
            # important to treat as Dictionary, since lmfit params & parinfo both have .items
            parnames,parvals = zip(*pars.items())
            parnames = [p.lower() for p in parnames]
            parvals = [p.value for p in parvals]
        elif parnames is None:
            parvals = pars
            parnames = self.parnames
        else:
            parvals = pars
        if len(pars) != len(parnames):
            # this should only be needed when other codes are changing the number of peaks
            # during a copy, as opposed to letting them be set by a __call__
            # (n_modelfuncs = n_ammonia can be called directly)
            # n_modelfuncs doesn't care how many peaks there are
            if len(pars) % len(parnames) == 0:
                parnames = [p for ii in range(len(pars)//len(parnames)) for p in parnames]
                npeaks = int(len(parvals) / npars)
                log.debug("Setting npeaks={0} npars={1}".format(npeaks, npars))
            else:
                raise ValueError("Wrong array lengths passed to n_ammonia!")

        self._components = []
        def L(x):
            v = np.zeros(len(x))
            for jj in range(int(npeaks)):
                modelkwargs = kwargs.copy()
                for ii in range(int(npars)):
                    name = parnames[ii+jj*int(npars)].strip('0123456789').lower()
                    modelkwargs.update({name:parvals[ii+jj*int(npars)]})
                v += self.modelfunc(x,**modelkwargs)
            return v
        return L

    def components(self, xarr, pars, hyperfine=False,
                   return_hyperfine_components=False, **kwargs):
        """
        Ammonia components don't follow the default, since in Galactic
        astronomy the hyperfine components should be well-separated.
        If you want to see the individual components overlaid, you'll need to
        pass hyperfine to the plot_fit call
        """

        comps=[]
        for ii in range(self.npeaks):
            if hyperfine or return_hyperfine_components:
                modelkwargs = dict(zip(self.parnames[ii*self.npars:(ii+1)*self.npars],
                                       pars[ii*self.npars:(ii+1)*self.npars]))
                comps.append(self.modelfunc(xarr, return_components=True,
                                            **modelkwargs))
            else:
                modelkwargs = dict(zip(self.parnames[ii*self.npars:(ii+1)*self.npars],
                                       pars[ii*self.npars:(ii+1)*self.npars]))
                comps.append([self.modelfunc(xarr, return_components=False,
                                             **modelkwargs)])

        modelcomponents = np.concatenate(comps)

        return modelcomponents


    def multinh3fit(self, xax, data, err=None,
                    parinfo=None,
                    quiet=True, shh=True,
                    debug=False,
                    maxiter=200,
                    use_lmfit=False,
                    veryverbose=False, **kwargs):
        """
        Fit multiple nh3 profiles (multiple can be 1)

        Parameters
        ----------
        xax : array
            x axis
        data : array
            y axis
        npeaks : int
            How many nh3 profiles to fit?  Default 1 (this could supersede onedgaussfit)
        err : array
            error corresponding to data
        params : list
            Fit parameters: [trot, tex, ntot (or tau), width, offset, ortho fraction] * npeaks
            If len(params) % 6 == 0, npeaks will be set to len(params) / 6.
            These parameters (and the related fixed, limited, min/max, names
            below) need to have length = 6*npeaks.  If npeaks > 1 and length =
            6, they will be replicated npeaks times, otherwise they will be
            reset to defaults:
        fixed : list
            Is parameter fixed?
        limitedmin : list 
        minpars : list
            set lower limits on each parameter (default: width>0, Tex and trot > Tcmb)
        limitedmax : list
        maxpars : list
            set upper limits on each parameter
        parnames : list
            default parameter names, important for setting kwargs in model
            ['trot','tex','ntot','width','xoff_v','fortho']
        quiet : bool
            should MPFIT output each iteration?
        shh : bool
            output final parameters?

        Returns
        -------
        mpp : model parameter object
            Fit parameters
        model : array
            The model array
        errors : array
            the fit errors
        chi2 : float
            the chi^2 value of the fit
        """

        if parinfo is None:
            parinfo = self.parinfo = self.make_parinfo(**kwargs)
        else:
            if isinstance(parinfo, ParinfoList):
                if not quiet:
                    log.info("Using user-specified parinfo.")
                self.parinfo = parinfo
            else:
                if not quiet:
                    log.info("Using something like a user-specified parinfo, but not.")
                self.parinfo = ParinfoList([p if isinstance(p,Parinfo) else Parinfo(p)
                                            for p in parinfo],
                                           preserve_order=True)

        fitfun_kwargs = dict((x,y) for (x,y) in kwargs.items()
                             if x not in ('npeaks', 'params', 'parnames',
                                          'fixed', 'limitedmin', 'limitedmax',
                                          'minpars', 'maxpars', 'tied',
                                          'max_tem_step'))
        fitfun_kwargs.update(self.modelfunc_kwargs)

        if 'use_lmfit' in fitfun_kwargs:
            raise KeyError("use_lmfit was specified in a location where it "
                           "is unacceptable")

        # not used: npars = len(parinfo)/self.npeaks

        self._validate_parinfo()

        def mpfitfun(x,y,err):
            if err is None:
                def f(p,fjac=None):
                    return [0,(y-self.n_ammonia(pars=p,
                                                parnames=parinfo.parnames,
                                                **fitfun_kwargs)(x))]
            else:
                def f(p,fjac=None):
                    return [0,(y-self.n_ammonia(pars=p,
                                                parnames=parinfo.parnames,
                                                **fitfun_kwargs)(x))/err]
            return f

        if veryverbose:
            log.info("GUESSES: ")
            log.info(str(parinfo))
            #log.info "\n".join(["%s: %s" % (p['parname'],p['value']) for p in parinfo])

        if use_lmfit:
            return self.lmfitter(xax, data, err=err,
                                 parinfo=parinfo,
                                 quiet=quiet,
                                 debug=debug)
        else:
            mp = mpfit(mpfitfun(xax,data,err),
                       parinfo=parinfo,
                       maxiter=maxiter,
                       quiet=quiet,
                       debug=debug)
            mpp = mp.params
            if mp.perror is not None:
                mpperr = mp.perror
            else:
                mpperr = mpp*0
            chi2 = mp.fnorm

        if mp.status == 0:
            raise Exception(mp.errmsg)

        for i,p in enumerate(mpp):
            parinfo[i]['value'] = p
            parinfo[i]['error'] = mpperr[i]

        if not shh:
            log.info("Fit status: {0}".format(mp.status))
            log.info("Fit message: {0}".format(mpfit_messages[mp.status]))
            log.info("Fit error message: {0}".format(mp.errmsg))
            log.info("Final fit values: ")
            for i,p in enumerate(mpp):
                log.info(" ".join((parinfo[i]['parname'], str(p), " +/- ",
                                   str(mpperr[i]))))
            log.info(" ".join(("Chi2: ", str(mp.fnorm)," Reduced Chi2: ",
                               str(mp.fnorm/len(data)), " DOF:",
                               str(len(data)-len(mpp)))))

        self.mp = mp

        self.parinfo = parinfo

        self.mpp = self.parinfo.values
        self.mpperr = self.parinfo.errors
        self.mppnames = self.parinfo.names
        self.model = self.n_ammonia(pars=self.mpp, parnames=self.mppnames,
                                    **fitfun_kwargs)(xax)

        indiv_parinfo = [self.parinfo[jj*self.npars:(jj+1)*self.npars]
                         for jj in range(int(len(self.parinfo)/self.npars))]
        modelkwargs = [dict([(p['parname'].strip("0123456789").lower(),
                              p['value']) for p in pi])
                       for pi in indiv_parinfo]
        self.tau_list = [self.modelfunc(xax, return_tau=True,**mk)
                         for mk in modelkwargs]

        return self.mpp,self.model,self.mpperr,chi2

    def moments(self, Xax, data, negamp=None, veryverbose=False,  **kwargs):
        """
        Returns a very simple and likely incorrect guess
        """

        # trot, TEX, ntot, width, center, ortho fraction
        return [20,10, 15, 1.0, 0.0, 1.0]

    def annotations(self):
        from decimal import Decimal # for formatting
        tex_key = {'trot':'T_R', 'tkin': 'T_K', 'tex':'T_{ex}', 'ntot':'N',
                   'fortho':'F_o', 'width':'\\sigma', 'xoff_v':'v',
                   'fillingfraction':'FF', 'tau':'\\tau_{1-1}',
                   'background_tb':'T_{BG}', 'delta':'T_R-T_{ex}'}
        # small hack below: don't quantize if error > value.  We want to see the values.
        label_list = []
        for pinfo in self.parinfo:
            parname = tex_key[pinfo['parname'].strip("0123456789").lower()]
            parnum = int(pinfo['parname'][-1])
            if pinfo['fixed']:
                formatted_value = "%s" % pinfo['value']
                pm = ""
                formatted_error=""
            else:
                formatted_value = Decimal("%g" % pinfo['value']).quantize(Decimal("%0.2g" % (min(pinfo['error'],pinfo['value']))))
                pm = "$\\pm$"
                formatted_error = Decimal("%g" % pinfo['error']).quantize(Decimal("%0.2g" % pinfo['error']))
            label = "$%s(%i)$=%8s %s %8s" % (parname, parnum, formatted_value, pm, formatted_error)
            label_list.append(label)
        labels = tuple(mpcb.flatten(label_list))
        return labels

    def make_parinfo(self, quiet=True,
                     npeaks=1,
                     params=(20,20,14,1.0,0.0,0.5), parnames=None,
                     fixed=(False,False,False,False,False,False),
                     limitedmin=(True,True,True,True,False,True),
                     limitedmax=(False,False,True,False,False,True),
                     minpars=(TCMB,TCMB,5,0,0,0),
                     maxpars=(0,0,25,0,0,1),
                     tied=('',)*6,
                     max_tem_step=1.,
                     **kwargs
                     ):

        if not quiet:
            log.info("Creating a 'parinfo' from guesses.")
        self.npars = int(len(params) / npeaks)

        if len(params) != npeaks and (len(params) / self.npars) > npeaks:
            npeaks = len(params) / self.npars
        npeaks = self.npeaks = int(npeaks)

        if isinstance(params,np.ndarray):
            params=params.tolist()
        # this is actually a hack, even though it's decently elegant
        # somehow, parnames was being changed WITHOUT being passed as a variable
        # this doesn't make sense - at all - but it happened.
        # (it is possible for self.parnames to have npars*npeaks elements where
        # npeaks > 1 coming into this function even though only 6 pars are specified;
        # _default_parnames is the workaround)
        if parnames is None:
            parnames = copy.copy(self._default_parnames)

        partype_dict = dict(zip(['params', 'parnames', 'fixed',
                                 'limitedmin', 'limitedmax', 'minpars',
                                 'maxpars', 'tied'],
                                [params, parnames, fixed, limitedmin,
                                 limitedmax, minpars, maxpars, tied]))

        # make sure all various things are the right length; if they're
        # not, fix them using the defaults
        # (you can put in guesses of length 12 but leave the rest length 6;
        # this code then doubles the length of everything else)
        for partype,parlist in iteritems(partype_dict):
            if len(parlist) != self.npars*self.npeaks:
                # if you leave the defaults, or enter something that can be
                # multiplied by npars to get to the right number of
                # gaussians, it will just replicate
                if len(parlist) == self.npars:
                    partype_dict[partype] *= npeaks
                elif len(parlist) > self.npars:
                    # DANGER:  THIS SHOULD NOT HAPPEN!
                    log.warning("WARNING!  Input parameters were longer than allowed for variable {0}".format(parlist))
                    partype_dict[partype] = partype_dict[partype][:self.npars]
                elif parlist==params: # this instance shouldn't really be possible
                    partype_dict[partype] = [20,20,1e10,1.0,0.0,0.5] * npeaks
                elif parlist==fixed:
                    partype_dict[partype] = [False] * len(params)
                elif parlist==limitedmax: # only fortho, fillingfraction have upper limits
                    partype_dict[partype] = (np.array(parnames) == 'fortho') + (np.array(parnames) == 'fillingfraction')
                elif parlist==limitedmin: # no physical values can be negative except velocity
                    partype_dict[partype] = (np.array(parnames) != 'xoff_v')
                elif parlist==minpars:
                    # all have minima of zero except kinetic temperature, which can't be below CMB.
                    # Excitation temperature technically can be, but not in this model
                    partype_dict[partype] = ((np.array(parnames) == 'trot') + (np.array(parnames) == 'tex')) * TCMB
                elif parlist==maxpars: # fractions have upper limits of 1.0
                    partype_dict[partype] = ((np.array(parnames) == 'fortho') + (np.array(parnames) == 'fillingfraction')).astype('float')
                elif parlist==parnames: # assumes the right number of parnames (essential)
                    partype_dict[partype] = list(parnames) * self.npeaks
                elif parlist==tied:
                    partype_dict[partype] = [_increment_string_number(t, ii*self.npars)
                                             for t in tied
                                             for ii in range(self.npeaks)]

        if len(parnames) != len(partype_dict['params']):
            raise ValueError("Wrong array lengths AFTER fixing them")

        # used in components.  Is this just a hack?
        self.parnames = partype_dict['parnames']

        parinfo = [{'n':ii, 'value':partype_dict['params'][ii],
                    'limits':[partype_dict['minpars'][ii],partype_dict['maxpars'][ii]],
                    'limited':[partype_dict['limitedmin'][ii],partype_dict['limitedmax'][ii]], 'fixed':partype_dict['fixed'][ii],
                    'parname':partype_dict['parnames'][ii]+str(int(ii/int(self.npars))),
                    'tied':partype_dict['tied'][ii],
                    'mpmaxstep':max_tem_step*float(partype_dict['parnames'][ii] in ('tex','trot')), # must force small steps in temperature (True = 1.0)
                    'error': 0}
                   for ii in range(len(partype_dict['params']))
                  ]

        # hack: remove 'fixed' pars
        #parinfo_with_fixed = parinfo
        #parinfo = [p for p in parinfo_with_fixed if not p['fixed']]
        #fixed_kwargs = dict((p['parname'].strip("0123456789").lower(),
        #                     p['value'])
        #                    for p in parinfo_with_fixed if p['fixed'])
        # don't do this - it breaks the NEXT call because npars != len(parnames) self.parnames = [p['parname'] for p in parinfo]
        # this is OK - not a permanent change
        #parnames = [p['parname'] for p in parinfo]
        # not OK self.npars = len(parinfo)/self.npeaks
        parinfo = ParinfoList([Parinfo(p) for p in parinfo], preserve_order=True)
        #import pdb; pdb.set_trace()
        return parinfo

    def _validate_parinfo(self,
                          must_be_limited={'trot': [True,False],
                                           'tex': [False,False],
                                           'ntot': [True, True],
                                           'width': [True, False],
                                           'xoff_v': [False, False],
                                           'tau': [False, False],
                                           'fortho': [True, True]},
                          required_limits={'trot': [0, None],
                                           'ntot': [5, 25],
                                           'width': [0, None],
                                           'fortho': [0,1]}):
        """
        Make sure the input parameters are all legitimate
        """
        for par in self.parinfo:
            limited = par.limited
            parname = par.parname.strip(string.digits).lower()
            mbl = must_be_limited[parname]

            for a,b,ul in zip(limited, mbl, ('a lower','an upper')):
                if b and not a:
                    raise ValueError("Parameter {0} must have {1} limit "
                                     "but no such limit is set.".format(
                                         parname, ul))
            if parname in required_limits:
                limits = par.limits
                rlimits = required_limits[parname]
                for a,b,op,ul in zip(limits, rlimits, (operator.lt,
                                                       operator.gt),
                                     ('a lower','an upper')):
                    if b is not None and op(a,b):
                        raise ValueError("Parameter {0} must have {1} limit "
                                         "at least {2} but it is set to {3}."
                                         .format(parname, ul, b, a))


    def parse_3par_guesses(self, guesses):
        """
        Try to convert a set of interactive guesses (peak, center, width) into
        guesses appropriate to the model.

        For NH3 models, we add in several extra parameters:
            tex = 2.73 * peak
            trot = tex * 2
            fortho = 0.5
            ntot = 15

        ntot is set to a constant ~10^15 because this results in optical depths
        near 1, so it forces the emission to be approximately significant.
        trot > tex so that we're in a physical regime to begin with.

        We assume tex = peak + 2.73 because most spectra are shown
        background-subtracted (single dish are always that way, interferometric
        data are intrinsically that way...) and otherwise the guessing will
        crash if you guess a number < 2.73.
        """
        gauss_npars = 3
        if len(guesses) % gauss_npars != 0:
            raise ValueError("Guesses passed to parse_3par_guesses must have "
                             "length % 3 == 0")
        npeaks = len(guesses) // gauss_npars
        npars = 6

        new_guesses = [-1, -1, 15, -1, -1, 0.5] * npeaks

        for ii in range(npeaks):
            peak = guesses[ii * gauss_npars + 0]
            center = guesses[ii * gauss_npars + 1]
            width = guesses[ii * gauss_npars + 2]
            new_guesses[ii*npars + 0] = (2.73+peak) * 2
            new_guesses[ii*npars + 1] = (2.73+peak)
            new_guesses[ii*npars + 3] = width
            new_guesses[ii*npars + 4] = center

        return new_guesses

class ammonia_model_vtau(ammonia_model):
    def __init__(self,
                 parnames=['trot', 'tex', 'tau', 'width', 'xoff_v', 'fortho'],
                 **kwargs):
        super(ammonia_model_vtau, self).__init__(parnames=parnames,
                                                 **kwargs)

    def moments(self, Xax, data, negamp=None, veryverbose=False,  **kwargs):
        """
        Returns a very simple and likely incorrect guess
        """

        # trot, TEX, ntot, width, center, ortho fraction
        return [20, 10, 10, 1.0, 0.0, 1.0]

    def _validate_parinfo(self,
                          must_be_limited={'trot': [True,False],
                                           'tex': [False,False],
                                           'tau': [True, False],
                                           'width': [True, False],
                                           'xoff_v': [False, False],
                                           'fortho': [True, True]},
                          required_limits={'trot': [0, None],
                                           'tex': [None,None],
                                           'width': [0, None],
                                           'tau': [0, None],
                                           'xoff_v': [None,None],
                                           'fortho': [0,1]}):
        supes = super(ammonia_model_vtau, self)
        supes._validate_parinfo(must_be_limited=must_be_limited,
                                required_limits=required_limits)
        return supes

    def make_parinfo(self,
                     params=(20,14,0.5,1.0,0.0,0.5),
                     fixed=(False,False,False,False,False,False),
                     limitedmin=(True,True,True,True,False,True),
                     limitedmax=(False,False,False,False,False,True),
                     minpars=(TCMB,TCMB,0,0,0,0),
                     maxpars=(0,0,0,0,0,1),
                     tied=('',)*6,
                     **kwargs
                     ):
        """
        parnames=['trot', 'tex', 'tau', 'width', 'xoff_v', 'fortho']
        """
        return super(ammonia_model_vtau, self).make_parinfo(params=params,
                                                            fixed=fixed,
                                                            limitedmax=limitedmax,
                                                            limitedmin=limitedmin,
                                                            minpars=minpars,
                                                            maxpars=maxpars,
                                                            tied=tied,
                                                            **kwargs)



class ammonia_model_vtau_thin(ammonia_model_vtau):
    def __init__(self,parnames=['tkin', 'tau', 'width', 'xoff_v', 'fortho'],
                 **kwargs):
        super(ammonia_model_vtau_thin, self).__init__(parnames=parnames,
                                                      npars=5,
                                                      **kwargs)
        self.modelfunc = ammonia_thin

    def _validate_parinfo(self,
                          must_be_limited={'tkin': [True,False],
                                           'tex': [False,False],
                                           'ntot': [True, True],
                                           'width': [True, False],
                                           'xoff_v': [False, False],
                                           'tau': [False, False],
                                           'fortho': [True, True]},
                          required_limits={'tkin': [0, None],
                                           'ntot': [5, 25],
                                           'width': [0, None],
                                           'fortho': [0,1]}):
        supes = super(ammonia_model_vtau_thin, self)
        return supes._validate_parinfo(must_be_limited=must_be_limited,
                                       required_limits=required_limits)


    def moments(self, Xax, data, negamp=None, veryverbose=False,  **kwargs):
        """
        Returns a very simple and likely incorrect guess
        """

        # trot, tau, width, center, ortho fraction
        return [20, 1, 1.0, 0.0, 1.0]

    def __call__(self,*args,**kwargs):
        return self.multinh3fit(*args, **kwargs)

    def make_parinfo(self,
                     params=(20,14,1.0,0.0,0.5),
                     fixed=(False,False,False,False,False),
                     limitedmin=(True,True,True,False,True),
                     limitedmax=(False,False,False,False,True),
                     minpars=(TCMB,0,0,0,0),
                     maxpars=(0,0,0,0,1),
                     tied=('',)*5,
                     **kwargs
                     ):
        return super(ammonia_model_vtau_thin, self).make_parinfo(params=params,
                                                                 fixed=fixed,
                                                                 limitedmax=limitedmax,
                                                                 limitedmin=limitedmin,
                                                                 minpars=minpars,
                                                                 maxpars=maxpars,
                                                                 tied=tied,
                                                                 **kwargs)

class ammonia_model_background(ammonia_model):
    def __init__(self,**kwargs):
        super(ammonia_model_background,self).__init__(npars=7,
                                                      parnames=['trot', 'tex',
                                                                'ntot',
                                                                'width',
                                                                'xoff_v',
                                                                'fortho',
                                                                'background_tb'])

    def moments(self, Xax, data, negamp=None, veryverbose=False,  **kwargs):
        """
        Returns a very simple and likely incorrect guess
        """

        # trot, TEX, ntot, width, center, ortho fraction
        return [20,10, 10, 1.0, 0.0, 1.0, TCMB]

    def make_parinfo(self, npeaks=1, err=None,
                     params=(20,20,14,1.0,0.0,0.5,TCMB), parnames=None,
                     fixed=(False,False,False,False,False,False,True),
                     limitedmin=(True,True,True,True,False,True,True),
                     limitedmax=(False,False,False,False,False,True,True),
                     minpars=(TCMB,TCMB,0,0,0,0,TCMB), parinfo=None,
                     maxpars=(0,0,0,0,0,1,TCMB),
                     tied=('',)*7,
                     quiet=True, shh=True,
                     veryverbose=False, **kwargs):
        return super(ammonia_model_background,
                     self).make_parinfo(npeaks=npeaks, err=err, params=params,
                                        parnames=parnames, fixed=fixed,
                                        limitedmin=limitedmin,
                                        limitedmax=limitedmax, minpars=minpars,
                                        parinfo=parinfo, maxpars=maxpars,
                                        tied=tied, quiet=quiet, shh=shh,
                                        veryverbose=veryverbose, **kwargs)

    def multinh3fit(self, xax, data, npeaks=1, err=None,
                    params=(20,20,14,1.0,0.0,0.5,TCMB), parnames=None,
                    fixed=(False,False,False,False,False,False,True),
                    limitedmin=(True,True,True,True,False,True,True),
                    limitedmax=(False,False,False,False,False,True,True),
                    minpars=(TCMB,TCMB,0,0,0,0,TCMB), parinfo=None,
                    maxpars=(0,0,0,0,0,1,TCMB),
                    tied=('',)*7,
                    quiet=True, shh=True,
                    veryverbose=False, **kwargs):
        return super(ammonia_model_background,
                     self).multinh3fit(xax, data, npeaks=npeaks, err=err,
                                       params=params, parnames=parnames,
                                       fixed=fixed, limitedmin=limitedmin,
                                       limitedmax=limitedmax, minpars=minpars,
                                       parinfo=parinfo, maxpars=maxpars,
                                       tied=tied, quiet=quiet, shh=shh,
                                       veryverbose=veryverbose, **kwargs)


class cold_ammonia_model(ammonia_model):
    def __init__(self,
                 parnames=['tkin', 'tex', 'ntot', 'width', 'xoff_v', 'fortho'],
                 **kwargs):
        super(cold_ammonia_model, self).__init__(parnames=parnames, **kwargs)
        self.modelfunc = cold_ammonia

    def _validate_parinfo(self,
                          must_be_limited={'tkin': [True,False],
                                           'tex': [False,False],
                                           'ntot': [True, False],
                                           'width': [True, False],
                                           'xoff_v': [False, False],
                                           'fortho': [True, True]},
                          required_limits={'tkin': [0, None],
                                           'tex': [None,None],
                                           'width': [0, None],
                                           'ntot': [0, None],
                                           'xoff_v': [None,None],
                                           'fortho': [0,1]}):
        supes = super(cold_ammonia_model, self)
        return supes._validate_parinfo(must_be_limited=must_be_limited,
                                       required_limits=required_limits)

class ammonia_model_restricted_tex(ammonia_model):
    """
    Ammonia model with an explicitly restricted excitation temperature
    such that tex <= trot, set by the "delta" parameter (tex = trot - delta)
    with delta > 0.  You can choose the ammonia funciton when you initialize
    it (e.g., ``ammonia_model_restricted_tex(ammonia_func=ammonia)`` or
    ``ammonia_model_restricted_tex(ammonia_func=cold_ammonia)``)
    """
    def __init__(self,
                 parnames=['trot', 'tex', 'ntot', 'width', 'xoff_v', 'fortho',
                           'delta'],
                 ammonia_func=ammonia,
                 **kwargs):
        super(ammonia_model_restricted_tex, self).__init__(npars=7,
                                                           parnames=parnames,
                                                           **kwargs)
        def ammonia_dtex(*args, **kwargs):
            """
            Strip out the 'delta' keyword
            """
            # for py2 compatibility, must strip out manually
            delta = kwargs.pop('delta') if 'delta' in kwargs else None
            try:
                np.testing.assert_allclose(kwargs['trot'] - kwargs['tex'],
                                        delta)
            except AssertionError:
                raise ValueError('trot - tex != delta, even though that was specified.'
                                 f"trot={kwargs['trot']}, tex={kwargs['tex']}, delta={delta}, trot-text={kwargs['trot']-kwargs['tex']}")
            return ammonia_func(*args, **kwargs)
        self.modelfunc = ammonia_dtex

    def n_ammonia(self, pars=None, parnames=None, **kwargs):
        if parnames is not None:
            for ii,pn in enumerate(parnames):
                if ii % 7 == 1 and 'tex' not in pn:
                    raise ValueError('bad parameter names')
                if ii % 7 == 6 and 'delta' not in pn:
                    raise ValueError('bad parameter names')
        if pars is not None:
            assert len(pars) % 7 == 0
            for ii in range(int(len(pars)/7)):
                try:
                    # Case A: they're param objects
                    # (setting the param directly can result in recursion errors)
                    pars[1+ii*7].value = pars[0+ii*7].value - pars[6+ii*7].value
                except AttributeError:
                    # Case B: they're just lists of values
                    pars[1+ii*7] = pars[0+ii*7] - pars[6+ii*7]

        supes = super(ammonia_model_restricted_tex, self)
        return supes.n_ammonia(pars=pars, parnames=parnames, **kwargs)


    def _validate_parinfo(self,
                          must_be_limited={'trot': [True,False],
                                           'tex': [False,False],
                                           'ntot': [True, False],
                                           'width': [True, False],
                                           'xoff_v': [False, False],
                                           'fortho': [True, True],
                                           'delta': [True, False],
                                          },
                          required_limits={'trot': [0, None],
                                           'tex': [None,None],
                                           'width': [0, None],
                                           'ntot': [0, None],
                                           'xoff_v': [None,None],
                                           'fortho': [0,1],
                                           'delta': [0, None],
                                          }):
        supes = super(ammonia_model_restricted_tex, self)
        return supes._validate_parinfo(must_be_limited=must_be_limited,
                                       required_limits=required_limits)

    def make_parinfo(self,
                     params=(20,20,0.5,1.0,0.0,0.5,0),
                     fixed=(False,False,False,False,False,False,False),
                     limitedmin=(True,True,True,True,False,True,True),
                     limitedmax=(False,False,False,False,False,True,False),
                     minpars=(TCMB,TCMB,0,0,0,0,0),
                     maxpars=(0,0,0,0,0,1,0),
                     tied=('','p[0]-p[6]','','','','',''),
                     **kwargs
                     ):
        """
        parnames=['trot', 'tex', 'ntot', 'width', 'xoff_v', 'fortho', 'delta']

        'delta' is the difference between tex and trot
        """
        supes = super(ammonia_model_restricted_tex, self)
        return supes.make_parinfo(params=params, fixed=fixed,
                                  limitedmax=limitedmax, limitedmin=limitedmin,
                                  minpars=minpars, maxpars=maxpars, tied=tied,
                                  **kwargs)






def _increment_string_number(st, count):
    """
    Increment a number in a string

    Expects input of the form: p[6]
    """

    import re
    dig = re.compile('[0-9]+')

    if dig.search(st):
        n = int(dig.search(st).group())
        result = dig.sub(str(n+count), st)
        return result
    else:
        return st
