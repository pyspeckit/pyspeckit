"""
========================================
Ammonia inversion transition TKIN fitter
========================================

Ammonia inversion transition TKIN fitter translated from Erik Rosolowsky's
http://svn.ok.ubc.ca/svn/signals/nh3fit/

.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>

Module API
^^^^^^^^^^

"""
import numpy as np
from ...mpfit import mpfit
from ...spectrum.parinfo import ParinfoList,Parinfo
from . import fitter
from . import model
import matplotlib.cbook as mpcb
import copy
from astropy import log
import astropy.units as u
from astropy import constants
from astropy.extern.six import iteritems
from . import mpfit_messages
import operator
import string

from .ammonia_constants import (line_names, freq_dict, aval_dict, ortho_dict,
                                voff_lines_dict, tau_wts_dict)

def ammonia(xarr, tkin=20, tex=None, ntot=14, width=1, xoff_v=0.0,
            fortho=0.0, tau=None, fillingfraction=None, return_tau=False,
            background_tb=2.7315,
            thin=False, verbose=False, return_components=False, debug=False,
            line_names=line_names):
    """
    Generate a model Ammonia spectrum based on input temperatures, column, and
    gaussian parameters

    Parameters
    ----------
    xarr: `pyspeckit.spectrum.units.SpectroscopicAxis`
        Array of wavelength/frequency values
    tex: float or None
        Excitation temperature. Assumed LTE if unspecified (``None``), if
        tex>tkin, or if ``thin`` is specified.
    ntot: float
        Total log column density of NH3.  Can be specified as a float in the
        range 5-25
    width: float
        Line width in km/s
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
    thin: bool
        uses a different parametetrization and requires only the optical depth,
        width, offset, and tkin to be specified.  In the 'thin' approximation,
        tex is not used in computation of the partition function - LTE is
        implicitly assumed
    return_components: bool
        Return a list of arrays, one for each hyperfine component, instead of
        just one array
    background_tb : float
        The background brightness temperature.  Defaults to TCMB.
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

    # Convert X-units to frequency in GHz
    xarr = xarr.as_unit('GHz')

    if tex is not None:
        # Yes, you certainly can have nonthermal excitation, tex>tkin.
        #if tex > tkin: # cannot have Tex > Tkin
        #    tex = tkin
        if thin: # tex is not used in this case
            tex = tkin
    else:
        tex = tkin

    if thin:
        ntot = 15
        if tau is None:
            raise ValueError("When using the 'thin' approximation, tau must "
                             "be used as a parameter.")
    elif 5 < ntot < 25:
        # allow ntot to be specified as a logarithm.  This is
        # safe because ntot < 1e10 gives a spectrum of all zeros, and the
        # plausible range of columns is not outside the specified range
        ntot = 10**ntot
    else:
        raise ValueError("ntot, the logarithmic total column density,"
                         " must be in the range 5 - 25")

    # fillingfraction is an arbitrary scaling for the data
    # The model will be (normal model) * fillingfraction
    if fillingfraction is None:
        fillingfraction = 1.0

    #ckms = 2.99792458e5
    ckms = constants.c.to(u.km/u.s).value
    ccms = constants.c.to(u.cm/u.s).value
    #Degeneracies
    # g1 = 1
    # g2 = 1
    #h = 6.6260693e-27
    h = constants.h.cgs.value
    #kb = 1.3806505e-16
    kb = constants.k_B.cgs.value
    # Dipole Moment in cgs (1.476 Debeye)
    #mu0 = 1.476e-18

    # Generate Partition Functions
    nlevs = 51
    jv=np.arange(nlevs)
    ortho = jv % 3 == 0
    para = ~ortho
    Jpara = jv[para]
    Jortho = jv[ortho]
    Brot = 298117.06e6
    Crot = 186726.36e6

    runspec = np.zeros(len(xarr))

    tau_dict = {}
    para_count = 0
    ortho_count = 1 # ignore 0-0

    if tau is not None and thin:
        """
        Use optical depth in the 1-1 line as a free parameter
        The optical depths of the other lines are then set by the kinetic temperature
        Tex is still a free parameter in the final spectrum calculation at the bottom
        (technically, I think this process assumes LTE; Tex should come into play in
        these equations, not just the final one)
        """
        dT0 = 41.5                    # Energy diff between (2,2) and (1,1) in K
        trot = tkin/(1+tkin/dT0*np.log(1+0.6*np.exp(-15.7/tkin)))
        tau_dict['oneone']     = tau
        tau_dict['twotwo']     = tau*(23.722/23.694)**2*4/3.*5/3.*np.exp(-41.5/trot)
        tau_dict['threethree'] = tau*(23.8701279/23.694)**2*3/2.*14./3.*np.exp(-101.1/trot)
        tau_dict['fourfour']   = tau*(24.1394169/23.694)**2*8/5.*9/3.*np.exp(-177.34/trot)
        line_names = tau_dict.keys()
    else:
        """
        Column density is the free parameter.  It is used in conjunction with
        the full partition function to compute the optical depth in each band
        Given the complexity of these equations, it would be worth my while to
        comment each step carefully.
        """
        Zpara = (2*Jpara+1)*np.exp(-h*(Brot*Jpara*(Jpara+1)+
                                       (Crot-Brot)*Jpara**2)/(kb*tkin))
        Zortho = 2*(2*Jortho+1)*np.exp(-h*(Brot*Jortho*(Jortho+1)+
                                           (Crot-Brot)*Jortho**2)/(kb*tkin))
        for linename in line_names:
            if ortho_dict[linename]:
                orthoparafrac = fortho
                Z = Zortho
                count = ortho_count
                ortho_count += 1
            else:
                orthoparafrac = 1.0-fortho
                Z = Zpara
                count = para_count # need to treat partition function separately
                para_count += 1

            # short variable names for readability
            frq = freq_dict[linename]
            partition = Z[count]
            aval = aval_dict[linename]

            # Friesen 2009 eqn A4 points out that the partition function actually says
            # how many molecules are in the NH3(1-1) state, both upper *and* lower.
            # population_upperlower = ntot * orthoparafrac * partition/(Z.sum())
            # population_upperstate = population_upperlower / (1+np.exp(h*frq/(kb*tex)))
            #
            # Note Jan 1, 2015: This is accounted for in the eqn below.  The
            # only difference is that I have used Tkin where Friesen et al 2009
            # use Tex.  Since Tex describes which states are populated, that may
            # be the correct one to use.

            # Total population of the higher energy inversion transition
            population_upperstate = ntot * orthoparafrac * partition/(Z.sum())

            tau_dict[linename] = (population_upperstate /
                                  (1. + np.exp(-h*frq/(kb*tkin) ))*ccms**2 /
                                  (8*np.pi*frq**2) * aval *
                                  (1-np.exp(-h*frq/(kb*tex))) /
                                  (width/ckms*frq*np.sqrt(2*np.pi)) )

    # allow tau(11) to be specified instead of ntot
    # in the thin case, this is not needed: ntot plays no role
    # this process allows you to specify tau without using the approximate equations specified
    # above.  It should remove ntot from the calculations anyway...
    if tau is not None and not thin:
        tau11_temp = tau_dict['oneone']
        # re-scale all optical depths so that tau is as specified, but the relative taus
        # are sest by the kinetic temperature and partition functions
        for linename,t in iteritems(tau_dict):
            tau_dict[linename] = t * tau/tau11_temp

    components =[]
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
            tauprof += (tau_dict[linename] * tau_wts[kk] *
                        np.exp(-(xarr.value+nuo-lines[kk])**2 /
                               (2.0*nuwidth[kk]**2)) * fillingfraction)
            components.append( tauprof )

        T0 = (h*xarr.value*1e9/kb) # "temperature" of wavelength
        if tau is not None and thin:
            #runspec = tauprof+runspec
            # is there ever a case where you want to ignore the optical depth function? I think no
            runspec = (T0/(np.exp(T0/tex)-1)-T0/(np.exp(T0/background_tb)-1))*(1-np.exp(-tauprof))+runspec
        else:
            runspec = (T0/(np.exp(T0/tex)-1)-T0/(np.exp(T0/background_tb)-1))*(1-np.exp(-tauprof))+runspec
        if runspec.min() < 0 and background_tb == 2.7315:
            raise ValueError("Model dropped below zero.  That is not possible normally.  Here are the input values: "+
                    ("tex: %f " % tex) +
                    ("tkin: %f " % tkin) +
                    ("ntot: %f " % ntot) +
                    ("width: %f " % width) +
                    ("xoff_v: %f " % xoff_v) +
                    ("fortho: %f " % fortho)
                    )

    if verbose or debug:
        log.info("tkin: %g  tex: %g  ntot: %g  width: %g  xoff_v: %g  fortho: %g  fillingfraction: %g" % (tkin,tex,ntot,width,xoff_v,fortho,fillingfraction))

    if return_components:
        return (T0/(np.exp(T0/tex)-1)-T0/(np.exp(T0/background_tb)-1))*(1-np.exp(-1*np.array(components)))

    if return_tau:
        return tau_dict

    return runspec


class ammonia_model(model.SpectralModel):

    def __init__(self,npeaks=1,npars=6,
                 parnames=['tkin','tex','ntot','width','xoff_v','fortho'],
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
                par.limits = (max(par.limits[0],2.73), par.limits[1])
            if 'tkin' in par.parname.lower():
                par.limited = (True,par.limited[1])
                par.limits = (max(par.limits[0],2.73), par.limits[1])
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
        self.fitunits = 'GHz'

    def __call__(self,*args,**kwargs):
        return self.multinh3fit(*args,**kwargs)

    def n_ammonia(self, pars=None, parnames=None, **kwargs):
        """
        Returns a function that sums over N ammonia line profiles, where N is the length of
        tkin,tex,ntot,width,xoff_v,fortho *OR* N = len(pars) / 6

        The background "height" is assumed to be zero (you must "baseline" your
        spectrum before fitting)

        *pars* [ list ]
            a list with len(pars) = (6-nfixed)n, assuming
            tkin,tex,ntot,width,xoff_v,fortho repeated

        *parnames* [ list ]
            len(parnames) must = len(pars).  parnames determine how the ammonia
            function parses the arguments
        """
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
                parnames = [p for ii in range(len(pars)/len(parnames)) for p in parnames]
                npars = len(parvals) / self.npeaks
            else:
                raise ValueError("Wrong array lengths passed to n_ammonia!")
        else:
            npars = int(len(parvals) / self.npeaks)

        self._components = []
        def L(x):
            v = np.zeros(len(x))
            for jj in range(int(self.npeaks)):
                modelkwargs = kwargs.copy()
                for ii in range(int(npars)):
                    name = parnames[ii+jj*int(npars)].strip('0123456789').lower()
                    modelkwargs.update({name:parvals[ii+jj*int(npars)]})
                v += ammonia(x,**modelkwargs)
            return v
        return L

    def components(self, xarr, pars, hyperfine=False, **kwargs):
        """
        Ammonia components don't follow the default, since in Galactic astronomy the hyperfine components should be well-separated.
        If you want to see the individual components overlaid, you'll need to pass hyperfine to the plot_fit call
        """

        comps=[]
        for ii in range(self.npeaks):
            if hyperfine:
                modelkwargs = dict(zip(self.parnames[ii*self.npars:(ii+1)*self.npars],pars[ii*self.npars:(ii+1)*self.npars]))
                comps.append( ammonia(xarr,return_components=True,**modelkwargs) )
            else:
                modelkwargs = dict(zip(self.parnames[ii*self.npars:(ii+1)*self.npars],pars[ii*self.npars:(ii+1)*self.npars]))
                comps.append( [ammonia(xarr,return_components=False,**modelkwargs)] )

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

        Inputs:
           xax - x axis
           data - y axis
           npeaks - How many nh3 profiles to fit?  Default 1 (this could supersede onedgaussfit)
           err - error corresponding to data

         These parameters need to have length = 6*npeaks.  If npeaks > 1 and length = 6, they will
         be replicated npeaks times, otherwise they will be reset to defaults:
           params - Fit parameters: [tkin, tex, ntot (or tau), width, offset, ortho fraction] * npeaks
                  If len(params) % 6 == 0, npeaks will be set to len(params) / 6
           fixed - Is parameter fixed?
           limitedmin/minpars - set lower limits on each parameter (default: width>0, Tex and Tkin > Tcmb)
           limitedmax/maxpars - set upper limits on each parameter
           parnames - default parameter names, important for setting kwargs in model ['tkin','tex','ntot','width','xoff_v','fortho']

           quiet - should MPFIT output each iteration?
           shh - output final parameters?

        Returns:
           Fit parameters
           Model
           Fit errors
           chi2
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
        if 'use_lmfit' in fitfun_kwargs:
            raise KeyError("use_lmfit was specified in a location where it "
                           "is unacceptable")

        npars = len(parinfo)/self.npeaks

        self._validate_parinfo()

        def mpfitfun(x,y,err):
            if err is None:
                def f(p,fjac=None): return [0,(y-self.n_ammonia(pars=p,
                                                                parnames=parinfo.parnames,
                                                                **fitfun_kwargs)(x))]
            else:
                def f(p,fjac=None): return [0,(y-self.n_ammonia(pars=p,
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
            if mp.perror is not None: mpperr = mp.perror
            else: mpperr = mpp*0
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
        modelkwargs = [
                dict([(p['parname'].strip("0123456789").lower(),p['value']) for p in pi])
                for pi in indiv_parinfo]
        self.tau_list = [ammonia(xax,return_tau=True,**mk) for mk in modelkwargs]

        return self.mpp,self.model,self.mpperr,chi2

    def moments(self, Xax, data, negamp=None, veryverbose=False,  **kwargs):
        """
        Returns a very simple and likely incorrect guess
        """

        # TKIN, TEX, ntot, width, center, ortho fraction
        return [20,10, 15, 1.0, 0.0, 1.0]

    def annotations(self):
        from decimal import Decimal # for formatting
        tex_key = {'tkin':'T_K', 'tex':'T_{ex}', 'ntot':'N', 'fortho':'F_o',
                   'width':'\\sigma', 'xoff_v':'v', 'fillingfraction':'FF',
                   'tau':'\\tau_{1-1}'}
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
                formatted_value = "%g" % (Decimal("%g" % pinfo['value']).quantize(Decimal("%0.2g" % (min(pinfo['error'],pinfo['value'])))))
                pm = "$\\pm$"
                formatted_error = "%g" % (Decimal("%g" % pinfo['error']).quantize(Decimal("%0.2g" % pinfo['error'])))
            label = "$%s(%i)$=%8s %s %8s" % (parname, parnum, formatted_value, pm, formatted_error)
            label_list.append(label)
        labels = tuple(mpcb.flatten(label_list))
        return labels

    def make_parinfo(self, quiet=True,
                     npeaks=1,
                     params=(20,20,14,1.0,0.0,0.5), parnames=None,
                     fixed=(False,False,False,False,False,False),
                     limitedmin=(True,True,True,True,False,True),
                     limitedmax=(False,False,False,False,False,True),
                     minpars=(2.73,2.73,0,0,0,0),
                     maxpars=(0,0,0,0,0,1),
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
        if parnames is None: parnames = copy.copy(self._default_parnames)

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
                    log.warn("WARNING!  Input parameters were longer than allowed for variable {0}".format(parlist))
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
                    partype_dict[partype] = ((np.array(parnames) == 'tkin') + (np.array(parnames) == 'tex')) * 2.73
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

        parinfo = [ {'n':ii, 'value':partype_dict['params'][ii],
                     'limits':[partype_dict['minpars'][ii],partype_dict['maxpars'][ii]],
                     'limited':[partype_dict['limitedmin'][ii],partype_dict['limitedmax'][ii]], 'fixed':partype_dict['fixed'][ii],
                     'parname':partype_dict['parnames'][ii]+str(int(ii/int(self.npars))),
                     'tied':partype_dict['tied'][ii],
                     'mpmaxstep':max_tem_step*float(partype_dict['parnames'][ii] in ('tex','tkin')), # must force small steps in temperature (True = 1.0)
                     'error': 0}
            for ii in range(len(partype_dict['params'])) ]

        # hack: remove 'fixed' pars
        #parinfo_with_fixed = parinfo
        #parinfo = [p for p in parinfo_with_fixed if not p['fixed']]
        #fixed_kwargs = dict((p['parname'].strip("0123456789").lower(),
        #                     p['value'])
        #                    for p in parinfo_with_fixed if p['fixed'])
        ## don't do this - it breaks the NEXT call because npars != len(parnames) self.parnames = [p['parname'] for p in parinfo]
        ## this is OK - not a permanent change
        #parnames = [p['parname'] for p in parinfo]
        ## not OK self.npars = len(parinfo)/self.npeaks
        parinfo = ParinfoList([Parinfo(p) for p in parinfo], preserve_order=True)
        #import pdb; pdb.set_trace()
        return parinfo

    def _validate_parinfo(self):
        """
        Make sure the input parameters are all legitimate
        """
        must_be_limited = {'tkin': [True,False],
                           'tex': [False,False],
                           'ntot': [True, False],
                           'width': [True, False],
                           'xoff_v': [False, False],
                           'tau': [False, False],
                           'fortho': [True, True]}
        required_limits = {'tkin': [0, None],
                           'ntot': [0, None],
                           'width': [0, None],
                           'fortho': [0,1]}
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

class ammonia_model_vtau(ammonia_model):
    def __init__(self,**kwargs):
        super(ammonia_model_vtau,self).__init__(parnames=['tkin','tex','tau','width','xoff_v','fortho'])

    def moments(self, Xax, data, negamp=None, veryverbose=False,  **kwargs):
        """
        Returns a very simple and likely incorrect guess
        """

        # TKIN, TEX, ntot, width, center, ortho fraction
        return [20, 10, 10, 1.0, 0.0, 1.0]

    def __call__(self,*args,**kwargs):
        return self.multinh3fit(*args,**kwargs)


class ammonia_model_vtau_thin(ammonia_model):
    def __init__(self,**kwargs):
        super(ammonia_model_vtau_thin, self).__init__(parnames=['tkin', 'tau',
                                                                'width',
                                                                'xoff_v',
                                                                'fortho'],
                                                      npars=5)

    def moments(self, Xax, data, negamp=None, veryverbose=False,  **kwargs):
        """
        Returns a very simple and likely incorrect guess
        """

        # TKIN, tau, width, center, ortho fraction
        return [20, 1, 1.0, 0.0, 1.0]

    def __call__(self,*args,**kwargs):
        return self.multinh3fit(*args, thin=True, **kwargs)

    def make_parinfo(self,
                     params=(20,14,1.0,0.0,0.5),
                     fixed=(False,False,False,False,False),
                     limitedmin=(True,True,True,False,True),
                     limitedmax=(False,False,False,False,True),
                     minpars=(2.73,0,0,0,0),
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
                                                      parnames=['tkin', 'tex',
                                                                'ntot',
                                                                'width',
                                                                'xoff_v',
                                                                'fortho',
                                                                'background_tb'])

    def moments(self, Xax, data, negamp=None, veryverbose=False,  **kwargs):
        """
        Returns a very simple and likely incorrect guess
        """

        # TKIN, TEX, ntot, width, center, ortho fraction
        return [20,10, 10, 1.0, 0.0, 1.0, 2.73]

    def __call__(self,*args,**kwargs):
        #if self.multisingle == 'single':
        #    return self.onepeakammoniafit(*args,**kwargs)
        #elif self.multisingle == 'multi':
        #    # Why is tied 6 instead of 7?
        return self.multinh3fit(*args,**kwargs)

    def make_parinfo(self, npeaks=1, err=None,
                    params=(20,20,14,1.0,0.0,0.5,2.73), parnames=None,
                    fixed=(False,False,False,False,False,False,True),
                    limitedmin=(True,True,True,True,False,True,True),
                    limitedmax=(False,False,False,False,False,True,True),
                    minpars=(2.73,2.73,0,0,0,0,2.73), parinfo=None,
                    maxpars=(0,0,0,0,0,1,2.73),
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
                    params=(20,20,14,1.0,0.0,0.5,2.73), parnames=None,
                    fixed=(False,False,False,False,False,False,True),
                    limitedmin=(True,True,True,True,False,True,True),
                    limitedmax=(False,False,False,False,False,True,True),
                    minpars=(2.73,2.73,0,0,0,0,2.73), parinfo=None,
                    maxpars=(0,0,0,0,0,1,2.73),
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

    def annotations(self):
        from decimal import Decimal # for formatting
        tex_key = {'tkin':'T_K', 'tex':'T_{ex}', 'ntot':'N', 'fortho':'F_o',
                   'width':'\\sigma', 'xoff_v':'v', 'fillingfraction':'FF',
                   'tau':'\\tau_{1-1}', 'background_tb':'T_{BG}'}
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
            label =  "$%s(%i)$=%8s %s %8s" % (parname, parnum, formatted_value, pm, formatted_error)
            label_list.append(label)
        labels = tuple(mpcb.flatten(label_list))
        return labels

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
