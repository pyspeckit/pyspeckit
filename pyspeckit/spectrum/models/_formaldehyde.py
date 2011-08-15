"""
Ammonia inversion transition TKIN fitter translated from Erik Rosolowsky's
http://svn.ok.ubc.ca/svn/signals/nh3fit/
"""
import numpy as np
from mpfit import mpfit
from .. import units
import matplotlib.cbook as mpcb

line_names = ['oneone','twotwo','threethree']
line_names = ['oneone_f10','oneone_f01','oneone_f22','oneone_f21','oneone_f12','oneone_f11']

# http://articles.adsabs.harvard.edu/abs/1971ApJ...169..429T has the most accurate freqs
freq_dict = { 
    'oneone':     4.82965996e9,
    'twotwo':     14.48848e9,
    'threethree': 28.97480e9,
    }
relative_strength_theory={
        'oneone_f10': 4,
        'oneone_f01': 4,
        'oneone_f22':15,
        'oneone_f21': 5,
        'oneone_f12': 5,
        'oneone_f11': 3,
        'twotwo_f11':1,
        'twotwo_f12':1,
        'twotwo_f21':1,
        'twotwo_f32':1,
        'twotwo_f33':1,
        'twotwo_f22':1,
        'twotwo_f23':1,
        'threethree_f22':1,
        'threethree_f44':1,
        'threethree_f33':1,
        }
hf_freq_dict={
        'oneone_f10':4.82965996e9 - 18.53e3,
        'oneone_f01':4.82965996e9 - 1.34e3,
        'oneone_f22':4.82965996e9 - 0.35e3,
        'oneone_f21':4.82965996e9 + 4.05e3,
        'oneone_f12':4.82965996e9 + 6.48e3,
        'oneone_f11':4.82965996e9 + 11.08e3,
        'twotwo_f11':14.48846e9,
        'twotwo_f12':14.48847e9,
        'twotwo_f21':14.48848e9,
        'twotwo_f32':14.48848e9,
        'twotwo_f33':14.48848e9,
        'twotwo_f22':14.48849e9,
        'twotwo_f23':14.48849e9,
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
voff_lines_dict={
        'oneone_f10': 18.53e3/4.82965996e9 * units.speedoflight_ms / 1000.0,
        'oneone_f01': 1.34e3 /4.82965996e9 * units.speedoflight_ms / 1000.0,
        'oneone_f22': 0.35e3 /4.82965996e9 * units.speedoflight_ms / 1000.0,
        'oneone_f21': 4.05e3 /4.82965996e9 * units.speedoflight_ms / 1000.0,
        'oneone_f12': 6.48e3 /4.82965996e9 * units.speedoflight_ms / 1000.0,
        'oneone_f11': 11.08e3/4.82965996e9 * units.speedoflight_ms / 1000.0,
        'twotwo_f11':14.48846e9,
        'twotwo_f12':14.48847e9,
        'twotwo_f21':14.48848e9,
        'twotwo_f32':14.48848e9,
        'twotwo_f33':14.48848e9,
        'twotwo_f22':14.48849e9,
        'twotwo_f23':14.48849e9,
        'threethree_f22':28.97478e9,
        'threethree_f44':28.97480e9,
        'threethree_f33':28.97481e9,
        }



class formaldehyde_model(object):

    def __init__(self):
        self.npeaks = 1
        self.npars = 6
        pass

    def formaldehyde(self, xarr, xunits='GHz', amp=1.0, width=1.0,
            xoff_v=0.0, line='oneone'):
        """
        Generate a model Formaldehyde spectrum based on input temperatures, column, and
        gaussian parameters


        (not implemented) if tau11 is specified, Ntot is ignored
        """

        # Convert X-units to frequency in GHz
        if xunits in units.frequency_dict:
            xarr = np.copy(xarr) * units.frequency_dict[xunits] / units.frequency_dict['GHz']
        elif xunits in units.velocity_dict:
            if line in freq_dict:
                xarr = (freq_dict[line] - (np.copy(xarr) * 
                        (units.velocity_dict[xunits] / units.velocity_dict['m/s'] / units.speedoflight_ms) *
                        freq_dict[line]) ) / units.frequency_dict['GHz']
            else:
                raise Exception("Xunits is velocity-type (%s) but line %s is not in the list." % (xunits,line))
        else:
            raise Exception("xunits not recognized: %s" % (xunits))

        ckms = 2.99792458e5
        ccms = ckms*1e5
        g1 = 1                
        g2 = 1                
        h = 6.6260693e-27     
        kb = 1.3806505e-16     

        runspec = np.zeros(len(xarr))

        for linename in line_names:
            voff_lines = np.array(voff_lines_dict[linename])
      
            lines = (1-voff_lines/ckms)*freq_dict[linename]
            nuwidth = np.abs(width/ckms*lines)
            nuoff = xoff_v/ckms*lines
      
            # strength array
            runspec += (1-relative_strength_theory[linename]*amp*\
                    np.exp(-(xarr+nuoff-freq_dict[linename])**2/(2*nuwidth**2)))
      
        return runspec

    def n_formaldehyde(self, pars=None, fittau=False, **kwargs):
        """
        Returns a function that sums over N ammonia line profiles, where N is the length of
        tkin,tex,Ntot,width,xoff_v,fortho *OR* N = len(pars) / 6

        The background "height" is assumed to be zero (you must "baseline" your
        spectrum before fitting)

        pars  - a list with len(pars) = 6n, assuming tkin,tex,Ntot,width,xoff_v,fortho repeated
        """
        if len(pars) % 6 == 0:
            tkin = [pars[ii] for ii in xrange(0,len(pars),6)]
            tex = [pars[ii] for ii in xrange(1,len(pars),6)]
            Ntot = [pars[ii] for ii in xrange(2,len(pars),6)]
            width = [pars[ii] for ii in xrange(3,len(pars),6)]
            xoff_v = [pars[ii] for ii in xrange(4,len(pars),6)]
            fortho = [pars[ii] for ii in xrange(5,len(pars),6)]
        elif not(len(tkin) == len(tex) == len(Ntot) == len(xoff_v) == len(width) == len(fortho)):
            raise ValueError("Wrong array lengths!")

        modelkwargs = kwargs.copy()
        def L(x):
            v = np.zeros(len(x))
            for i in range(len(tkin)):
                modelkwargs.update({'tkin':tkin[i], 'tex':tex[i],
                        'width':width[i], 'xoff_v':xoff_v[i],
                        'fortho':fortho[i]})
                if fittau:
                    modelkwargs.update({'tau11':Ntot[i]})
                else:
                    modelkwargs.update({'Ntot':Ntot[i]})
                v += self.ammonia(x,**modelkwargs)
            return v
        return L

    def multinh3fit(self, xax, data, npeaks=1, err=None, params=[20,20,1e10,1.0,0.0,0.5],
            fixed=[False,False,False,False,False,False],
            limitedmin=[True,True,True,True,False,True],
            limitedmax=[False,False,False,False,False,True], minpars=[2.73,2.73,0,0,0,0],
            maxpars=[0,0,0,0,0,1], quiet=True, shh=True, veryverbose=False, **kwargs):
        """
        Fit multiple nh3 profiles

        Inputs:
           xax - x axis
           data - y axis
           npeaks - How many nh3 profiles to fit?  Default 1 (this could supersede onedgaussfit)
           err - error corresponding to data

         These parameters need to have length = 6*npeaks.  If npeaks > 1 and length = 6, they will
         be replicated npeaks times, otherwise they will be reset to defaults:
           params - Fit parameters: [amplitude, offset, Gfwhm, Lfwhm] * npeaks
                  If len(params) % 6 == 0, npeaks will be set to len(params) / 6
           fixed - Is parameter fixed?
           limitedmin/minpars - set lower limits on each parameter (default: width>0, Tex and Tkin > Tcmb)
           limitedmax/maxpars - set upper limits on each parameter

           quiet - should MPFIT output each iteration?
           shh - output final parameters?

        Returns:
           Fit parameters
           Model
           Fit errors
           chi2
        """

        self.npars = 6

        if len(params) != npeaks and (len(params) / self.npars) > npeaks:
            npeaks = len(params) / self.npars 
        self.npeaks = npeaks

        if isinstance(params,np.ndarray): params=params.tolist()

        # make sure all various things are the right length; if they're not, fix them using the defaults
        for parlist in (params,fixed,limitedmin,limitedmax,minpars,maxpars):
            if len(parlist) != self.npars*self.npeaks:
                # if you leave the defaults, or enter something that can be multiplied by 3 to get to the
                # right number of gaussians, it will just replicate
                if len(parlist) == self.npars: 
                    parlist *= npeaks 
                elif parlist==params:
                    parlist[:] = [20,20,1e10,1.0,0.0,0.5] * npeaks
                elif parlist==fixed:
                    parlist[:] = [False,False,False,False,False,False] * npeaks
                elif parlist==limitedmax:
                    parlist[:] = [False,False,False,False,False,True] * npeaks
                elif parlist==limitedmin:
                    parlist[:] = [True,True,True,True,False,True] * npeaks
                elif parlist==minpars:
                    parlist[:] = [2.73,0,0,0,0,0] * npeaks
                elif parlist==maxpars:
                    parlist[:] = [0,0,0,0,0,1] * npeaks

        def mpfitfun(x,y,err):
            if err is None:
                def f(p,fjac=None): return [0,(y-self.n_ammonia(pars=p, **kwargs)(x))]
            else:
                def f(p,fjac=None): return [0,(y-self.n_ammonia(pars=p, **kwargs)(x))/err]
            return f

        parnames = {0:"TKIN",1:"TEX",2:"NTOT",3:"WIDTH",4:"XOFF_V",5:"FORTHO"}

        parinfo = [ {'n':ii, 'value':params[ii],
            'limits':[minpars[ii],maxpars[ii]],
            'limited':[limitedmin[ii],limitedmax[ii]], 'fixed':fixed[ii],
            'parname':parnames[ii%self.npars]+str(ii/self.npars), 
            'mpmaxstep':0,'error':ii} 
            for ii in xrange(len(params)) ]
        parinfo[0]['mpmaxstep'] = 1.0
        parinfo[1]['mpmaxstep'] = 1.0

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

        if not shh:
            print "Fit message: ",mp.errmsg
            print "Final fit values: "
            for i,p in enumerate(mpp):
                parinfo[i]['value'] = p
                print parinfo[i]['parname'],p," +/- ",mpperr[i]
            print "Chi2: ",mp.fnorm," Reduced Chi2: ",mp.fnorm/len(data)," DOF:",len(data)-len(mpp)

        if mpp[1] > mpp[0]: mpp[1] = mpp[0]  # force Tex>Tkin to Tex=Tkin (already done in n_ammonia)
        self.mp = mp
        self.mpp = mpp
        self.mpperr = mpperr
        self.model = self.n_ammonia(pars=mpp,**kwargs)(xax)
        return mpp,self.n_ammonia(pars=mpp,**kwargs)(xax),mpperr,chi2

    __call__ = multinh3fit

    def moments(self, Xax, data, negamp=None, veryverbose=False,  **kwargs):
        """
        Returns a very simple and likely incorrect guess
        """

        # TKIN, TEX, NTOT, width, center, ortho fraction
        return [20,10, 1e15, 1.0, 0.0, 1.0]

    def annotations(self):
        label_list = [ (
                "$T_K(%i)$=%6.4g $\\pm$ %6.4g" % (jj,self.mpp[0+jj*self.npars],self.mpperr[0+jj*self.npars]),
                "$T_{ex}(%i)$=%6.4g $\\pm$ %6.4g" % (jj,self.mpp[1+jj*self.npars],self.mpperr[1+jj*self.npars]),
                "$N$(%i)=%6.4g $\\pm$ %6.4g" % (jj,self.mpp[2+jj*self.npars],self.mpperr[2+jj*self.npars]),
                "$\\sigma(%i)$=%6.4g $\\pm$ %6.4g" % (jj,self.mpp[3+jj*self.npars],self.mpperr[3+jj*self.npars]),
                "$v(%i)$=%6.4g $\\pm$ %6.4g" % (jj,self.mpp[4+jj*self.npars],self.mpperr[4+jj*self.npars]),
                "$F_o(%i)$=%6.4g $\\pm$ %6.4g" % (jj,self.mpp[5+jj*self.npars],self.mpperr[5+jj*self.npars])
                          ) for jj in range(self.npeaks)]
        labels = tuple(mpcb.flatten(label_list))
        return labels
