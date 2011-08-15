"""
Ammonia inversion transition TKIN fitter translated from Erik Rosolowsky's
http://svn.ok.ubc.ca/svn/signals/nh3fit/
"""
import numpy as np
from mpfit import mpfit
from .. import units
from . import fitter,model
import matplotlib.cbook as mpcb
import copy

line_names = ['oneone','twotwo','threethree']
line_names = ['oneone_f10','oneone_f01','oneone_f22','oneone_f21','oneone_f12','oneone_f11','twotwo_f11','twotwo_f12','twotwo_f21','twotwo_f32','twotwo_f33','twotwo_f22','twotwo_f23']

# http://adsabs.harvard.edu/abs/1971ApJ...169..429T has the most accurate freqs
# http://adsabs.harvard.edu/abs/1972ApJ...174..463T [twotwo]
freq_dict = { 
    'oneone':     4.82965996e9,
    'twotwo':     14.48847881e9,
    'threethree': 28.97480e9,
    }
relative_strength_theory={
        'oneone_f10':  4.,
        'oneone_f01':  4.,
        'oneone_f22': 15.,
        'oneone_f21':  5.,
        'oneone_f12':  5.,
        'oneone_f11':  3.,
        'twotwo_f11': 15.,
        'twotwo_f12':  5.,
        'twotwo_f21':  5.,
        'twotwo_f32': 5.19,
        'twotwo_f33': 41.48,
        'twotwo_f22': 23.15,
        'twotwo_f23': 5.19,
        'threethree_f22':1,
        'threethree_f44':1,
        'threethree_f33':1,
        }
relative_strength_total_degeneracy={
        'oneone_f10': 36.,
        'oneone_f01': 36.,
        'oneone_f22': 36.,
        'oneone_f21': 36.,
        'oneone_f12': 36.,
        'oneone_f11': 36.,
        'twotwo_f11': 100.01,
        'twotwo_f12': 100.01,
        'twotwo_f21': 100.01,
        'twotwo_f32': 100.01,
        'twotwo_f33': 100.01,
        'twotwo_f22': 100.01,
        'twotwo_f23': 100.01,
        'threethree_f22':3.0,
        'threethree_f44':3.0,
        'threethree_f33':3.0,
        }
hf_freq_dict={
        'oneone_f10':4.82965996e9 - 18.53e3,
        'oneone_f01':4.82965996e9 - 1.34e3,
        'oneone_f22':4.82965996e9 - 0.35e3,
        'oneone_f21':4.82965996e9 + 4.05e3,
        'oneone_f12':4.82965996e9 + 6.48e3,
        'oneone_f11':4.82965996e9 + 11.08e3,
        'twotwo_f11':14.48847881e9 - 19.97e3,
        'twotwo_f12':14.48847881e9 -  7.03e3,
        'twotwo_f21':14.48847881e9 -  2.20e3,
        'twotwo_f32':14.48847881e9 +  0.12e3,
        'twotwo_f33':14.48847881e9 +  0.89e3,
        'twotwo_f22':14.48847881e9 + 10.74e3,
        'twotwo_f23':14.48847881e9 + 11.51e3,
        'threethree_f22':28.97478e9,
        'threethree_f44':28.97480e9,
        'threethree_f33':28.97481e9,
        }
freq_dict.update(hf_freq_dict)
aval_dict = {
    'oneone':     10**-8.44801,  #64*!pi**4/(3*h*c**3)*nu11**3*mu0**2*(1/2.)
    'twotwo':     10**-7.49373,  #64*!pi**4/(3*h*c**3)*nu22**3*mu0**2*(2/3.)
    'threethree': 10**-6.89179,  #64*!pi**4/(3*h*c**3)*nu33**3*mu0**2*(3/4.)
    }
hf_aval_dict={
        'oneone_f10':10**-8.92509,
        'oneone_f01':10**-8.44797,
        'oneone_f22':10**-8.57294,
        'oneone_f21':10**-9.05004,
        'oneone_f12':10**-8.82819,
        'oneone_f11':10**-9.05009,
        'twotwo_f11':10**-7.61876,
        'twotwo_f12':10**-8.09586,
        'twotwo_f21':10**-8.31771,
        'twotwo_f32':10**-8.44804,
        'twotwo_f33':10**-7.54494,
        'twotwo_f22':10**-7.65221,
        'twotwo_f23':10**-8.30191,
        'threethree_f22':10**-6.94294,
        'threethree_f44':10**-6.91981,
        'threethree_f33':10**-6.96736,
        }
ortho_dict = {
    'oneone':     False,
    'twotwo':     False,
    'threethree': False,
    }
n_ortho = np.arange(0,28,3) # 0..3..27
n_para = np.array([x for x in range(28) if x % 3 != 0])

voff_lines_dict = {
        'oneone': [(hf_freq_dict[f]-freq_dict['oneone'])/freq_dict['oneone']*units.speedoflight_ms for f in hf_freq_dict.keys() if "oneone" in f],
        'twotwo': [(hf_freq_dict[f]-freq_dict['twotwo'])/freq_dict['twotwo']*units.speedoflight_ms for f in hf_freq_dict.keys() if "twotwo" in f],
        'threethree': [(hf_freq_dict[f]-freq_dict['threethree'])/freq_dict['threethree']*units.speedoflight_ms for f in hf_freq_dict.keys() if "threethree" in f],
        }
voff_lines_dict={ # opposite signs of freq offset
        'oneone_f10': + 18.53e3/freq_dict['oneone'] * units.speedoflight_ms / 1000.0,
        'oneone_f01': + 1.34e3 /freq_dict['oneone'] * units.speedoflight_ms / 1000.0,
        'oneone_f22': + 0.35e3 /freq_dict['oneone'] * units.speedoflight_ms / 1000.0,
        'oneone_f21': - 4.05e3 /freq_dict['oneone'] * units.speedoflight_ms / 1000.0,
        'oneone_f12': - 6.48e3 /freq_dict['oneone'] * units.speedoflight_ms / 1000.0,
        'oneone_f11': - 11.08e3/freq_dict['oneone'] * units.speedoflight_ms / 1000.0,
        'twotwo_f11': + 19.97e3/freq_dict['twotwo'] * units.speedoflight_ms / 1000.0,
        'twotwo_f12': +  7.03e3/freq_dict['twotwo'] * units.speedoflight_ms / 1000.0,
        'twotwo_f21': +  2.20e3/freq_dict['twotwo'] * units.speedoflight_ms / 1000.0,
        'twotwo_f32': -  0.12e3/freq_dict['twotwo'] * units.speedoflight_ms / 1000.0,
        'twotwo_f33': -  0.89e3/freq_dict['twotwo'] * units.speedoflight_ms / 1000.0,
        'twotwo_f22': - 10.74e3/freq_dict['twotwo'] * units.speedoflight_ms / 1000.0,
        'twotwo_f23': - 11.51e3/freq_dict['twotwo'] * units.speedoflight_ms / 1000.0,
        'threethree_f22':28.97478e9,
        'threethree_f44':28.97480e9,
        'threethree_f33':28.97481e9,
        }


def formaldehyde_vtau(xarr, Tex=1.0, xoff_v=0.0, width=1.0, tau=1.0 , 
        return_components=False, Tbackground=2.73 ):
    """
    Generate a model Formaldehyde spectrum based on simple gaussian parameters with
    ampltiude set by Tex and tau. Fitted parameters will depend on the optical depth
    of hyperfine components

    """

    # Convert X-units to frequency in GHz
    xarr = copy.copy(xarr)
    xarr.convert_to_unit('Hz', quiet=True)

    ckms = 2.99792458e5

    tau_nu_cumul = np.zeros(len(xarr))
    if np.any(np.isnan((tau,Tex,width,xoff_v))):
        if return_components:
            return [tau_nu_cumul] * len(line_names)
        else:
            return tau_nu_cumul

    components =[]
    for linename in line_names:
        voff_lines = np.array(voff_lines_dict[linename])
  
        lines = (1-voff_lines/ckms)*freq_dict[linename]
        if width == 0:
            tau_nu = xarr*0
        else:
            nuwidth = np.abs(width/ckms*lines)
            nuoff = xoff_v/ckms*lines
            # the total optical depth, which is being fitted, should be the sum of the components
            tau_line = (tau * relative_strength_theory[linename]) / relative_strength_total_degeneracy[linename]
      
            tau_nu = np.array(tau_line * np.exp(-(xarr+nuoff-freq_dict[linename])**2/(2.0*nuwidth**2)))
            tau_nu[tau_nu!=tau_nu] = 0 # avoid nans
            components.append( tau_nu )
        tau_nu_cumul += tau_nu

    # add a list of the individual 'component' spectra to the total components...

    if return_components:
        return (1.0-np.exp(-np.array(components)))*(Tex-Tbackground)

    spec = (1.0-np.exp(-np.array(tau_nu_cumul)))*(Tex-Tbackground)
  
    return spec

formaldehyde_vtau_fitter = model.SpectralModel(formaldehyde_vtau,4,
        parnames=['Tex','center','width','tau'], 
        parlimited=[(False,False),(False,False), (True,False), (True,False)], 
        parlimits=[(0,0), (0,0), (0,0), (0,0)],
        fitunits='Hz' )

formaldehyde_vtau_vheight_fitter = model.SpectralModel(fitter.vheightmodel(formaldehyde_vtau),5,
        parnames=['height','Tex','center','width','tau'], 
        parlimited=[(False,False), (False,False),(False,False), (True,False), (True,False)], 
        parlimits=[(0,0), (0,0), (0,0), (0,0), (0,0)],
        fitunits='Hz' )


def formaldehyde(xarr, amp=1.0, xoff_v=0.0, width=1.0, 
        return_components=False ):
    """
    Generate a model Formaldehyde spectrum based on simple gaussian parameters

    """

    # Convert X-units to frequency in GHz
    xarr = copy.copy(xarr)
    xarr.convert_to_unit('Hz', quiet=True)

    ckms = 2.99792458e5

    runspec = np.zeros(len(xarr))
    if np.any(np.isnan((amp,width,xoff_v))):
        if return_components:
            return [runspec] * len(line_names)
        else:
            return runspec

    components =[]
    for linename in line_names:
        voff_lines = np.array(voff_lines_dict[linename])
  
        lines = (1-voff_lines/ckms)*freq_dict[linename]
        if width == 0:
            speccomp = xarr*0
        else:
            nuwidth = np.abs(width/ckms*lines)
            nuoff = xoff_v/ckms*lines
      
            speccomp = np.array(relative_strength_theory[linename] * np.exp(-(xarr+nuoff-freq_dict[linename])**2/(2.0*nuwidth**2)))
            speccomp[speccomp!=speccomp] = 0 # avoid nans
            components.append( speccomp )
        runspec += speccomp

    # add a list of the individual 'component' spectra to the total components...

    if return_components:
        return np.array(components)*amp/runspec.max()  

    runspec *= amp/runspec.max()
  
    return runspec

formaldehyde_fitter = model.SpectralModel(formaldehyde, 3,
        parnames=['amp','center','width'], 
        parlimited=[(False,False),(False,False), (True,False)], 
        parlimits=[(0,0), (0,0), (0,0)],
        fitunits='Hz' )

formaldehyde_vheight_fitter = model.SpectralModel(fitter.vheightmodel(formaldehyde), 4,
        parnames=['height','amp','center','width'], 
        parlimited=[(False,False),(False,False),(False,False), (True,False)], 
        parlimits=[(0,0), (0,0), (0,0), (0,0)],
        fitunits='Hz' )



class formaldehyde_model(fitter.SimpleFitter):

    def __init__(self,multisingle='multi'):
        self.npars = 3
        self.npeaks = 1

        self.onepeakformaldehydeian = fitter.vheightmodel(formaldehyde)
        self.onepeakformaldehydefit = self._fourparfitter(self.onepeakformaldehydeian)

        if multisingle in ('multi','single'):
            self.multisingle = multisingle
        else:
            raise Exception("multisingle must be multi or single")

    def __call__(self,*args,**kwargs):
        if self.multisingle == 'single':
            return self.onepeakformaldehydefit(*args,**kwargs)
        elif self.multisingle == 'multi':
            return self.multiformaldehydefit(*args,**kwargs)


    def n_formaldehyde(self, pars=None, modelfunc=formaldehyde, **kwargs):
        """
        Returns a function that sums over N h2co line profiles 
        The background "height" is assumed to be zero (you must "baseline" your
        spectrum before fitting)

        pars  - a list with len(pars) = 3n
        """
        if len(pars) % 3 == 0:
            a = [pars[ii] for ii in xrange(0,len(pars),3)]
            dx = [pars[ii] for ii in xrange(1,len(pars),3)]
            sigma = [pars[ii] for ii in xrange(2,len(pars),3)]
        elif not(len(tkin) == len(tex) == len(Ntot) == len(xoff_v) == len(width) == len(fortho)):
            raise ValueError("Wrong array lengths!")

        self._components = []
        modelkwargs = kwargs.copy()
        def L(x):
            v = np.zeros(len(x))
            for i in range(len(a)):
                modelkwargs.update({'amp':a[i], 'xoff_v':dx[i],
                        'width':sigma[i]})
                v += modelfunc(x,**modelkwargs)
            return v
        return L


    def multiformaldehydefit(self, xax, data, npeaks=1, err=None, params=[1,0,1],
            fixed=[False,False,False], limitedmin=[False,False,True],
            limitedmax=[False,False,False], minpars=[0,0,0], maxpars=[0,0,0],
            quiet=True, shh=True, veryverbose=False, negamp=None,
            tied = ['', '', ''], **kwargs):
        """
        An improvement on onepeakformaldehydefit.  Lets you fit multiple formaldehydeians.

        Inputs:
           xax - x axis [must be a SpectroscopicAxis instance]
           data - y axis
           npeaks - How many formaldehydeians to fit?  Default 1 (this could supersede onepeakformaldehydefit)
           err - error corresponding to data

         These parameters need to have length = 3*npeaks.  If npeaks > 1 and length = 3, they will
         be replicated npeaks times, otherwise they will be reset to defaults:
           params - Fit parameters: [amplitude, offset, width] * npeaks
                  If len(params) % 3 == 0, npeaks will be set to len(params) / 3
           fixed - Is parameter fixed?
           limitedmin/minpars - set lower limits on each parameter (default: width>0)
           limitedmax/maxpars - set upper limits on each parameter
           tied - link parameters together

           quiet - should MPFIT output each iteration?
           shh - output final parameters?

           kwargs are passed to mpfit

        Returns:
           Fit parameters
           Model
           Fit errors
           chi2
        """

        for p in params:
            if np.isnan(p):
                raise ValueError("ERROR: A parameter guess is NAN")

        if len(params) != npeaks and (len(params) / 3) > npeaks:
            self.npeaks = len(params) / 3 
        else:
            self.npeaks = npeaks

        if isinstance(params,np.ndarray): params=params.tolist()
        if not isinstance(xax,units.SpectroscopicAxis): 
            raise TypeError("X axis must be a SpectroscopicAxis instance.")

        # make sure all various things are the right length; if they're not, fix them using the defaults
        # multiformaldehydefit should process negamp directly if kwargs.has_key('negamp') is False: kwargs['negamp'] = None 
        pardict = {"params":params,"fixed":fixed,"limitedmin":limitedmin,"limitedmax":limitedmax,"minpars":minpars,"maxpars":maxpars,"tied":tied}
        for parlistname in pardict:
            parlist = pardict[parlistname]
            if len(parlist) != 3*self.npeaks:
                # if you leave the defaults, or enter something that can be multiplied by 3 to get to the
                # right number of formaldehydeians, it will just replicate
                if veryverbose: print "Correcting length of parameter %s" % parlistname
                if len(parlist) == 3: 
                    parlist *= self.npeaks 
                elif parlistname=="params":
                    parlist[:] = [1,0,1] * self.npeaks
                elif parlistname=="fixed":
                    parlist[:] = [False,False,False] * self.npeaks
                elif parlistname=="limitedmax":
                    if negamp is None: parlist[:] = [False,False,False] * self.npeaks
                    elif negamp is False: parlist[:] = [False,False,False] * self.npeaks
                    else: parlist[:] = [True,False,False] * self.npeaks
                elif parlistname=="limitedmin":
                    if negamp is None: parlist[:] = [False,False,True] * self.npeaks  # Lines can't have negative width!
                    elif negamp is False: parlist[:] = [True,False,True] * self.npeaks
                    else: parlist[:] = [False,False,True] * self.npeaks                   
                elif parlistname=="minpars" or parlistname=="maxpars":
                    parlist[:] = [0,0,0] * self.npeaks
                elif parlistname=="tied":
                    parlist[:] = ['','',''] * self.npeaks
                    
        # mpfit doesn't recognize negamp, so get rid of it now that we're done setting limitedmin/max and min/maxpars
        #if kwargs.has_key('negamp'): kwargs.pop('negamp')

        def mpfitfun(x,y,err):
            if err is None:
                def f(p,fjac=None): return [0,(y-self.n_formaldehyde(pars=p)(x))]
            else:
                def f(p,fjac=None): return [0,(y-self.n_formaldehyde(pars=p)(x))/err]
            return f

        parnames = {0:"AMPLITUDE",1:"SHIFT",2:"WIDTH"}

        parinfo = [ {'n':ii, 'value':params[ii],
            'limits':[minpars[ii],maxpars[ii]],
            'limited':[limitedmin[ii],limitedmax[ii]], 'fixed':fixed[ii],
            'parname':parnames[ii%3]+str(ii/3), 'error':ii, 'tied':tied[ii]} 
            for ii in xrange(len(params)) ]

        if veryverbose:
            print "GUESSES: "
            print "\n".join(["%s: %s" % (p['parname'],p['value']) for p in parinfo])

        mp = mpfit(mpfitfun(xax,data,err),parinfo=parinfo,quiet=quiet,**kwargs)
        mpp = mp.params
        if mp.perror is not None: mpperr = mp.perror
        else: mpperr = mpp*0
        chi2 = mp.fnorm

        if mp.status == 0:
            raise Exception(mp.errmsg)

        if not shh:
            print "Fit status: ",mp.status
            print "Fit error message: ",mp.errmsg
            print "Fit message: ",mpfit_messages[mp.status]
            print "Final fit values: "
            for i,p in enumerate(mpp):
                parinfo[i]['value'] = p
                print parinfo[i]['parname'],p," +/- ",mpperr[i]
            print "Chi2: ",mp.fnorm," Reduced Chi2: ",mp.fnorm/len(data)," DOF:",len(data)-len(mpp)

        self.mp = mp
        self.mpp = mpp
        self.mpperr = mpperr
        self.model = self.n_formaldehyde(pars=mpp)(xax)
        return mpp,self.n_formaldehyde(pars=mpp)(xax),mpperr,chi2

    def components(self, xarr, pars):

        modelcomponents = np.concatenate(
            [formaldehyde(xarr,
                amp = pars[3*i],
                xoff_v = pars[3*i+1],
                width = pars[3*i+2],
                return_components=True)
            for i in range(self.npeaks)])

        return modelcomponents

    def annotations(self):
        label_list = [(
                "$A(%i)$=%6.4g $\\pm$ %6.4g" % (jj,self.mpp[0+jj*self.npars],self.mpperr[0+jj*self.npars]),
                "$v(%i)$=%6.4g $\\pm$ %6.4g" % (jj,self.mpp[1+jj*self.npars],self.mpperr[1+jj*self.npars]),
                "$\\sigma(%i)$=%6.4g $\\pm$ %6.4g" % (jj,self.mpp[2+jj*self.npars],self.mpperr[2+jj*self.npars])
                          ) for jj in range(self.npeaks)]
        labels = tuple(mpcb.flatten(label_list))
        return labels

    def integral(self, modelpars, linename='oneone'):
        """
        Return the integral of the individual components (ignoring height)
        """
        # produced by directly computing the integral of gaussians and formaldehydeians as a function of 
        # line width and then fitting that with a broken logarithmic power law
        # The errors are <0.5% for all widths
        formaldehyde_to_gaussian_ratio_coefs = {
                'lt0.1_oneone': np.array([ -5.784020,-40.058798,-111.172706,-154.256411,-106.593122,-28.933119]),
                'gt0.1_oneone': np.array([  0.038548, -0.071162, -0.045710,  0.183828, -0.145429,  0.040039]),
                'lt0.1_twotwo': np.array([  1.156561,  6.638570, 11.782065, -0.429536,-24.860297,-27.902274, -9.510288]),
                'gt0.1_twotwo': np.array([ -0.090646,  0.078204,  0.123181, -0.175590,  0.089506, -0.034687,  0.008676]),
                }


        integ = 0
        if len(modelpars) % 3 == 0:
            for amp,cen,width in np.reshape(modelpars,[len(modelpars)/3,3]):
                gaussint = amp*width*np.sqrt(2.0*np.pi)
                cftype = "gt0.1_"+linename if width > 0.1 else "lt0.1_"+linename
                correction_factor = 10**np.polyval(formaldehyde_to_gaussian_ratio_coefs[cftype], np.log10(width) )
                # debug statement print "Two components of the integral: amp %g, width %g, gaussint %g, correction_factor %g " % (amp,width,gaussint,correction_factor)
                integ += gaussint*correction_factor

        return integ
