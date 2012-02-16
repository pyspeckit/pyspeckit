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
from pyspeckit.mpfit import mpfit
import fitter
import matplotlib.cbook as mpcb
import copy


line_names = ['oneone','twotwo','threethree','fourfour']

freq_dict = { 
    'oneone':     23.694506e9,
    'twotwo':     23.722633335e9,
    'threethree': 23.8701296e9,
    'fourfour':   24.1394169e9,
    }
aval_dict = {
    'oneone':     1.712e-7,  #64*!pi**4/(3*h*c**3)*nu11**3*mu0**2*(1/2.)
    'twotwo':     2.291e-7,  #64*!pi**4/(3*h*c**3)*nu22**3*mu0**2*(2/3.)
    'threethree': 2.625e-7,  #64*!pi**4/(3*h*c**3)*nu33**3*mu0**2*(3/4.)
    'fourfour':   3.167e-7,  #64*!pi**4/(3*h*c**3)*nu44**3*mu0**2*(4/5.)
    }
ortho_dict = {
    'oneone':     False,
    'twotwo':     False,
    'threethree': True,
    'fourfour':   False,
    }
n_ortho = np.arange(0,28,3) # 0..3..27
n_para = np.array([x for x in range(28) if x % 3 != 0])

voff_lines_dict = {
    'oneone': [19.8513, 19.3159, 7.88669, 7.46967, 7.35132, 0.460409, 0.322042,
        -0.0751680, -0.213003, 0.311034, 0.192266, -0.132382, -0.250923, -7.23349,
        -7.37280, -7.81526, -19.4117, -19.5500],
    'twotwo':[26.5263, 26.0111, 25.9505, 16.3917, 16.3793, 15.8642, 0.562503,
        0.528408, 0.523745, 0.0132820, -0.00379100, -0.0132820, -0.501831,
        -0.531340, -0.589080, -15.8547, -16.3698, -16.3822, -25.9505, -26.0111,
        -26.5263],
    'threethree':[29.195098, 29.044147, 28.941877, 28.911408, 21.234827,
        21.214619, 21.136387, 21.087456, 1.005122, 0.806082, 0.778062,
        0.628569, 0.016754, -0.005589, -0.013401, -0.639734, -0.744554,
        -1.031924, -21.125222, -21.203441, -21.223649, -21.076291, -28.908067,
        -28.938523, -29.040794, -29.191744],
    'fourfour':[  0.        , -30.49783692,  30.49783692,   0., 24.25907811,
        -24.25907811,   0.        ]
                      }

tau_wts_dict = {
    'oneone': [0.0740740, 0.148148, 0.0925930, 0.166667, 0.0185190, 0.0370370,
        0.0185190, 0.0185190, 0.0925930, 0.0333330, 0.300000, 0.466667,
        0.0333330, 0.0925930, 0.0185190, 0.166667, 0.0740740, 0.148148],
    'twotwo': [0.00418600, 0.0376740, 0.0209300, 0.0372090, 0.0260470,
        0.00186000, 0.0209300, 0.0116280, 0.0106310, 0.267442, 0.499668,
        0.146512, 0.0116280, 0.0106310, 0.0209300, 0.00186000, 0.0260470,
        0.0372090, 0.0209300, 0.0376740, 0.00418600],
    'threethree': [0.012263, 0.008409, 0.003434, 0.005494, 0.006652, 0.008852,
        0.004967, 0.011589, 0.019228, 0.010387, 0.010820, 0.009482, 0.293302,
        0.459109, 0.177372, 0.009482, 0.010820, 0.019228, 0.004967, 0.008852,
        0.006652, 0.011589, 0.005494, 0.003434, 0.008409, 0.012263],
    'fourfour': [0.2431, 0.0162, 0.0162, 0.3008, 0.0163, 0.0163, 0.3911]}

def ammonia(xarr, tkin=20, tex=None, ntot=1e14, width=1,
        xoff_v=0.0, fortho=0.0, tau=None, fillingfraction=None, return_tau=False,
        thin=False, verbose=False, return_components=False, debug=False ):
    """
    Generate a model Ammonia spectrum based on input temperatures, column, and
    gaussian parameters

    ntot can be specified as a column density (e.g., 10^15) or a log-column-density (e.g., 15)

    tex can be specified or can be assumed LTE if unspecified, if tex>tkin, or if "thin"
        is specified

    "thin" uses a different parametetrization and requires only the optical depth, width, offset,
        and tkin to be specified.  In the 'thin' approximation, tex is not used in computation of
        the partition function - LTE is implicitly assumed

    If tau is specified, ntot is NOT fit but is set to a fixed value
    fillingfraction is an arbitrary scaling factor to apply to the model
    fortho is the ortho/(ortho+para) fraction.  The default is to assume all ortho.
    xoff_v is the velocity offset in km/s 

    tau refers to the optical depth of the 1-1 line.  The optical depths of the
    other lines are fixed relative to tau_oneone

    (not implemented) if tau is specified, ntot is ignored
    """

    # Convert X-units to frequency in GHz
    xarr = xarr.as_unit('GHz')

    if tex is not None:
        if tex > tkin: # cannot have Tex > Tkin
            tex = tkin 
        elif thin: # tex is not used in this case
            tex = tkin
    else:
        tex = tkin

    if thin:
        ntot = 1e15
    elif 5 < ntot < 25: 
        # allow ntot to be specified as a logarithm.  This is
        # safe because ntot < 1e10 gives a spectrum of all zeros, and the
        # plausible range of columns is not outside the specified range
        ntot = 10**ntot
    elif (25 < ntot < 1e5) or (ntot < 5):
        # these are totally invalid for log/non-log
        return 0

    # fillingfraction is an arbitrary scaling for the data
    # The model will be (normal model) * fillingfraction
    if fillingfraction is None:
        fillingfraction = 1.0

    ckms = 2.99792458e5
    ccms = ckms*1e5
    g1 = 1                
    g2 = 1                
    h = 6.6260693e-27     
    kb = 1.3806505e-16     
    mu0 = 1.476e-18               # Dipole Moment in cgs (1.476 Debeye)
  
    # Generate Partition Functions  
    nlevs = 51
    jv=np.arange(nlevs)
    ortho = jv % 3 == 0
    para = True-ortho
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
            tau_dict[linename] = (ntot * orthoparafrac * Z[count]/(Z.sum()) / ( 1
                + np.exp(-h*freq_dict[linename]/(kb*tkin) )) * ccms**2 /
                (8*np.pi*freq_dict[linename]**2) * aval_dict[linename]*
                (1-np.exp(-h*freq_dict[linename]/(kb*tex))) /
                (width/ckms*freq_dict[linename]*np.sqrt(2*np.pi)) )

    # allow tau(11) to be specified instead of ntot
    # in the thin case, this is not needed: ntot plays no role
    # this process allows you to specify tau without using the approximate equations specified
    # above.  It should remove ntot from the calculations anyway...
    if tau is not None and not thin:
        tau11_temp = tau_dict['oneone']
        # re-scale all optical depths so that tau is as specified, but the relative taus
        # are sest by the kinetic temperature and partition functions
        for linename,t in tau_dict.iteritems():
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
        for kk,no in enumerate(nuoff):
            tauprof += (tau_dict[linename] * tau_wts[kk] *
                    np.exp(-(xarr+no-lines[kk])**2 / (2.0*nuwidth[kk]**2)) *
                    fillingfraction)
            components.append( tauprof )
  
        T0 = (h*xarr*1e9/kb) # "temperature" of wavelength
        if tau is not None and thin:
            #runspec = tauprof+runspec
            # is there ever a case where you want to ignore the optical depth function? I think no
            runspec = (T0/(np.exp(T0/tex)-1)-T0/(np.exp(T0/2.73)-1))*(1-np.exp(-tauprof))+runspec
        else:
            runspec = (T0/(np.exp(T0/tex)-1)-T0/(np.exp(T0/2.73)-1))*(1-np.exp(-tauprof))+runspec
        if runspec.min() < 0:
            raise ValueError("Model dropped below zero.  That is not possible normally.  Here are the input values: "+
                    ("tex: %f " % tex) + 
                    ("tkin: %f " % tkin) + 
                    ("ntot: %f " % ntot) + 
                    ("width: %f " % width) + 
                    ("xoff_v: %f " % xoff_v) + 
                    ("fortho: %f " % fortho)
                    )

    if verbose:
        print "tkin: %g  tex: %g  ntot: %g  width: %g  xoff_v: %g  fortho: %g  fillingfraction: %g" % (tkin,tex,ntot,width,xoff_v,fortho,fillingfraction)

    if return_components:
        return (T0/(np.exp(T0/tex)-1)-T0/(np.exp(T0/2.73)-1))*(1-np.exp(-1*np.array(components)))

    if return_tau:
        return tau_dict
  
    return runspec

class ammonia_model(fitter.SimpleFitter):

    def __init__(self,npeaks=1,npars=6,multisingle='multi',**kwargs):
        self.npeaks = npeaks
        self.npars = npars
        self._default_parnames = ['tkin','tex','ntot','width','xoff_v','fortho']
        self.parnames = copy.copy(self._default_parnames)

        # all fitters must have declared modelfuncs, which should take the fitted pars...
        self.modelfunc = ammonia
        self.n_modelfunc = self.n_ammonia

        self.onepeakammonia = fitter.vheightmodel(ammonia)
        #self.onepeakammoniafit = self._fourparfitter(self.onepeakammonia)

        if multisingle in ('multi','single'):
            self.multisingle = multisingle
        else:
            raise Exception("multisingle must be multi or single")

        self.modelfunc_kwargs = kwargs

    def __call__(self,*args,**kwargs):
        if self.multisingle == 'single':
            return self.onepeakammoniafit(*args,**kwargs)
        elif self.multisingle == 'multi':
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
        if parnames is None:
            parnames = self.parnames
        if len(pars) != len(parnames):
            raise ValueError("Wrong array lengths passed to n_ammonia!")
        else:
            npars = len(pars) / self.npeaks

        self._components = []
        def L(x):
            v = np.zeros(len(x))
            for jj in xrange(self.npeaks):
                modelkwargs = kwargs.copy()
                for ii in xrange(npars):
                    modelkwargs.update({parnames[ii+jj*npars].strip('0123456789'):pars[ii+jj*npars]})
                v += ammonia(x,**modelkwargs)
            return v
        return L

    def components(self, xarr, pars, hyperfine=False):
        """
        Ammonia components don't follow the default, since in Galactic astronomy the hyperfine components should be well-separated.
        If you want to see the individual components overlaid, you'll need to pass hyperfine to the plot_fit call
        """

        comps=[]
        for ii in xrange(self.npeaks):
            if hyperfine:
                modelkwargs = dict(zip(self.parnames[ii*self.npars:(ii+1)*self.npars],pars[ii*self.npars:(ii+1)*self.npars]))
                comps.append( ammonia(xarr,return_components=True,**modelkwargs) )
            else:
                modelkwargs = dict(zip(self.parnames[ii*self.npars:(ii+1)*self.npars],pars[ii*self.npars:(ii+1)*self.npars]))
                comps.append( [ammonia(xarr,return_components=False,**modelkwargs)] )

        modelcomponents = np.concatenate(comps)

        return modelcomponents


    def multinh3fit(self, xax, data, npeaks=1, err=None, 
            params=(20,20,14,1.0,0.0,0.5),
            parnames=None,
            fixed=(False,False,False,False,False,False),
            limitedmin=(True,True,True,True,False,True),
            limitedmax=(False,False,False,False,False,True), minpars=(2.73,2.73,0,0,0,0),
            parinfo=None,
            maxpars=(0,0,0,0,0,1), quiet=True, shh=True, veryverbose=False, **kwargs):
        """
        Fit multiple nh3 profiles

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
            self.npars = len(params) / npeaks

            if len(params) != npeaks and (len(params) / self.npars) > npeaks:
                npeaks = len(params) / self.npars 
            self.npeaks = npeaks

            if isinstance(params,np.ndarray): params=params.tolist()
            # this is actually a hack, even though it's decently elegant
            # somehow, parnames was being changed WITHOUT being passed as a variable
            # this doesn't make sense - at all - but it happened.
            # (it is possible for self.parnames to have npars*npeaks elements where
            # npeaks > 1 coming into this function even though only 6 pars are specified;
            # _default_parnames is the workaround)
            if parnames is None: parnames = copy.copy(self._default_parnames)

            partype_dict = dict(zip(['params','parnames','fixed','limitedmin','limitedmax','minpars','maxpars'],
                    [params,parnames,fixed,limitedmin,limitedmax,minpars,maxpars]))

            # make sure all various things are the right length; if they're not, fix them using the defaults
            for partype,parlist in partype_dict.iteritems():
                if len(parlist) != self.npars*self.npeaks:
                    # if you leave the defaults, or enter something that can be multiplied by npars to get to the
                    # right number of gaussians, it will just replicate
                    if len(parlist) == self.npars: 
                        partype_dict[partype] *= npeaks 
                    elif len(parlist) > self.npars:
                        # DANGER:  THIS SHOULD NOT HAPPEN!
                        print "WARNING!  Input parameters were longer than allowed for variable ",parlist
                        partype_dict[partype] = partype_dict[partype][:self.npars]
                    elif parlist==params: # this instance shouldn't really be possible
                        partype_dict[partype] = [20,20,1e10,1.0,0.0,0.5] * npeaks
                    elif parlist==fixed:
                        partype_dict[partype] = [False] * len(params)
                    elif parlist==limitedmax: # only fortho, fillingfraction have upper limits
                        partype_dict[partype] = (np.array(parnames) == 'fortho') + (np.array(parnames) == 'fillingfraction')
                    elif parlist==limitedmin: # no physical values can be negative except velocity
                        partype_dict[partype] = (np.array(parnames) != 'xoff_v')
                    elif parlist==minpars: # all have minima of zero except kinetic temperature, which can't be below CMB.  Excitation temperature technically can be, but not in this model
                        partype_dict[partype] = ((np.array(parnames) == 'tkin') + (np.array(parnames) == 'tex')) * 2.73
                    elif parlist==maxpars: # fractions have upper limits of 1.0
                        partype_dict[partype] = ((np.array(parnames) == 'fortho') + (np.array(parnames) == 'fillingfraction')).astype('float')
                    elif parlist==parnames: # assumes the right number of parnames (essential)
                        partype_dict[partype] = list(parnames) * self.npeaks 

            if len(parnames) != len(partype_dict['params']):
                raise ValueError("Wrong array lengths AFTER fixing them")

            # used in components.  Is this just a hack?
            self.parnames = partype_dict['parnames']

            parinfo = [ {'n':ii, 'value':partype_dict['params'][ii],
                'limits':[partype_dict['minpars'][ii],partype_dict['maxpars'][ii]],
                'limited':[partype_dict['limitedmin'][ii],partype_dict['limitedmax'][ii]], 'fixed':partype_dict['fixed'][ii],
                'parname':partype_dict['parnames'][ii]+str(ii/self.npars),
                'mpmaxstep':float(partype_dict['parnames'][ii] in ('tex','tkin')), # must force small steps in temperature (True = 1.0)
                'error': 0} 
                for ii in xrange(len(partype_dict['params'])) ]

            # hack: remove 'fixed' pars
            parinfo_with_fixed = parinfo
            parinfo = [p for p in parinfo_with_fixed if not p['fixed']]
            fixed_kwargs = dict((p['parname'].strip("0123456789"),p['value']) for p in parinfo_with_fixed if p['fixed'])
            # don't do this - it breaks the NEXT call because npars != len(parnames) self.parnames = [p['parname'] for p in parinfo]
            # this is OK - not a permanent change
            parnames = [p['parname'] for p in parinfo]
            # not OK self.npars = len(parinfo)/self.npeaks
        else: 
            self.parinfo = parinfo
            parinfo_with_fixed = None
            fixed_kwargs = {}

        fitfun_kwargs = dict(kwargs.items()+fixed_kwargs.items())

        npars = len(parinfo)/self.npeaks

        # (fortho0 is not fortho)
        # this doesn't work if parinfo_with_fixed is not None:
        # this doesn't work     for p in parinfo_with_fixed:
        # this doesn't work         # users can change the defaults while holding them fixed 
        # this doesn't work         if p['fixed']:
        # this doesn't work             kwargs.update({p['parname']:p['value']})

        def mpfitfun(x,y,err):
            if err is None:
                def f(p,fjac=None): return [0,(y-self.n_ammonia(pars=p, parnames=[pi['parname'] for pi in parinfo], **fitfun_kwargs)(x))]
            else:
                def f(p,fjac=None): return [0,(y-self.n_ammonia(pars=p, parnames=[pi['parname'] for pi in parinfo], **fitfun_kwargs)(x))/err]
            return f

        if veryverbose:
            print "GUESSES: "
            print "\n".join(["%s: %s" % (p['parname'],p['value']) for p in parinfo])

        mp = mpfit(mpfitfun(xax,data,err),parinfo=parinfo,quiet=quiet)
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
            print "Fit status: ",mp.status
            print "Fit message: ",mp.errmsg
            print "Final fit values: "
            for i,p in enumerate(mpp):
                print parinfo[i]['parname'],p," +/- ",mpperr[i]
            print "Chi2: ",mp.fnorm," Reduced Chi2: ",mp.fnorm/len(data)," DOF:",len(data)-len(mpp)

        if any(['tex' in s for s in parnames]) and any(['tkin' in s for s in parnames]):
            texnum = (i for i,s in enumerate(parnames) if 'tex' in s)
            tkinnum = (i for i,s in enumerate(parnames) if 'tkin' in s)
            for txn,tkn in zip(texnum,tkinnum):
                if mpp[txn] > mpp[tkn]: mpp[txn] = mpp[tkn]  # force Tex>Tkin to Tex=Tkin (already done in n_ammonia)
        self.mp = mp

        if parinfo_with_fixed is not None:
            # self self.parinfo preserving the 'fixed' parameters 
            self.parinfo = parinfo_with_fixed
            # ORDER MATTERS!
            for p in parinfo:
                self.parinfo[p['n']] = p
        else:
            self.parinfo = parinfo
        #self.parinfo = parinfo

        self.mpp = np.array([p['value'] for p in self.parinfo])
        self.mpperr = np.array([p['error'] for p in self.parinfo])
        self.model = self.n_ammonia(pars=self.mpp, parnames=[p['parname'] for p in self.parinfo], **kwargs)(xax)
        #if self.model.sum() == 0:
        #    print "DON'T FORGET TO REMOVE THIS ERROR!"
        #    raise ValueError("Model is zeros.")

        indiv_parinfo = [self.parinfo[jj*self.npars:(jj+1)*self.npars] for jj in xrange(len(self.parinfo)/self.npars)]
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
        return [20,10, 1e15, 1.0, 0.0, 1.0]

    def annotations(self):
        from decimal import Decimal # for formatting
        tex_key = {'tkin':'T_K','tex':'T_{ex}','ntot':'N','fortho':'F_o','width':'\\sigma','xoff_v':'v','fillingfraction':'FF','tau':'\\tau_{1-1}'}
        # small hack below: don't quantize if error > value.  We want to see the values.
        label_list = []
        for pinfo in self.parinfo:
            parname = tex_key[pinfo['parname'].strip("0123456789")]
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

class ammonia_model_vtau(ammonia_model):
    def __init__(self,**kwargs):
        super(ammonia_model_vtau,self).__init__()
        self.parnames = ['tkin','tex','tau','width','xoff_v','fortho']

    def moments(self, Xax, data, negamp=None, veryverbose=False,  **kwargs):
        """
        Returns a very simple and likely incorrect guess
        """

        # TKIN, TEX, ntot, width, center, ortho fraction
        return [20,10, 1, 1.0, 0.0, 1.0]

    def __call__(self,*args,**kwargs):
        if self.multisingle == 'single':
            return self.onepeakammoniafit(*args,**kwargs)
        elif self.multisingle == 'multi':
            return self.multinh3fit(*args,**kwargs)

