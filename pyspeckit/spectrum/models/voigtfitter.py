"""
====================
Voigt Profile Fitter
====================
"""
import numpy
from numpy.ma import median
from numpy import pi
from mpfit import mpfit
import matplotlib.cbook as mpcb
from . import fitter

try:
    import scipy.special
    scipyOK = True
except ImportError:
    scipyOK = False
    #print "Voigt profile fitting will fail.  Requires scipy.special"

class voigt_fitter(fitter.SimpleFitter):

    def __init__(self,multisingle='multi'):
        self.npars = 4
        self.npeaks = 1
        if multisingle in ('multi','single'):
            self.multisingle = multisingle
        else:
            raise Exception("multisingle must be multi or single")

    def __call__(self,*args,**kwargs):
        if self.multisingle == 'single':
            return self.onedvoigtfit(*args,**kwargs)
        elif self.multisingle == 'multi':
            return self.multivoigtfit(*args,**kwargs)

    def voigt(self,xarr,amp,xcen,Gfwhm,Lfwhm):
        """
        voigt profile

        V(x,sig,gam) = Re(w(z))/(sig*sqrt(2*pi))
        z = (x+i*gam)/(sig*sqrt(2))

        Converted from 
        http://mail.scipy.org/pipermail/scipy-user/2011-January/028327.html
        """

        if scipyOK:
            tmp = 1.0/scipy.special.wofz(numpy.zeros((len(xarr))) \
                  +1j*numpy.sqrt(numpy.log(2.0))*Lfwhm).real
            tmp = tmp*amp* \
                  scipy.special.wofz(2*numpy.sqrt(numpy.log(2.0))*(xarr-xcen)/Gfwhm+1j* \
                  numpy.sqrt(numpy.log(2.0))*Lfwhm).real
            return tmp
        else:
            raise ImportError("Couldn't import scipy, therefore cannot do voigt profile stuff")

    def n_voigt(self,pars=None):
        """
        Returns a function that sums over N voigt profiles, where N is the length of
        a,dx,Gfwhm,Lfwhm *OR* N = len(pars) / 4

        The background "height" is assumed to be zero (you must "baseline" your
        spectrum before fitting)

        pars  - a list with len(pars) = 3n, assuming a,dx,sigma repeated
        #dx    - offset (velocity center) values
        #width - line widths (Lorentzian FWHM)
        #a     - amplitudes
        """
        if len(pars) % 4 == 0:
            a = [pars[ii] for ii in xrange(0,len(pars),4)]
            dx = [pars[ii] for ii in xrange(1,len(pars),4)]
            Gfwhm = [pars[ii] for ii in xrange(2,len(pars),4)]
            Lfwhm = [pars[ii] for ii in xrange(3,len(pars),4)]
        elif not(len(dx) == len(width) == len(a)):
            raise ValueError("Wrong array lengths! dx: %i  width %i  a: %i" % (len(dx),len(width),len(a)))

        def L(x):
            v = numpy.zeros(len(x))
            for i in range(len(dx)):
                v += self.voigt(x,a[i],dx[i],Gfwhm[i],Lfwhm[i])
            return v
        return L

    def onedvoigtfit(self,xax, data, err=None,
            params=[0,1,0,1,1],fixed=[False,False,False,False,False],
            limitedmin=[False,False,False,True,True],
            limitedmax=[False,False,False,False,False], minpars=[0,0,0,0,0],
            maxpars=[0,0,0,0,0], quiet=True, shh=True, veryverbose=False,
            vheight=True, negamp=False, usemoments=False, **kwargs):
        """
        Inputs:
           xax - x axis
           data - y axis
           err - error corresponding to data

           params - Fit parameters: Height of background, Amplitude, Shift, Width
           fixed - Is parameter fixed?
           limitedmin/minpars - set lower limits on each parameter (default: width>0)
           limitedmax/maxpars - set upper limits on each parameter
           quiet - should MPFIT output each iteration?
           shh - output final parameters?
           usemoments - replace default parameters with moments

        Returns:
           Fit parameters
           Model
           Fit errors
           chi2
        """

        def mpfitfun(x,y,err):
            if err is None:
                def f(p,fjac=None): return [0,(y-self.n_voigt(pars=p[1:])(x)+p[0])]
            else:
                def f(p,fjac=None): return [0,(y-self.n_voigt(pars=p[1:])(x)+p[0])/err]
            return f

        if xax == None:
            xax = numpy.arange(len(data))

        if vheight is False: 
            height = params[0]
            fixed[0] = True
        if usemoments:
            params = moments(xax, data, vheight=vheight, negamp=negamp,
                    veryverbose=veryverbose)
            if vheight is False: params = [height]+params
            if veryverbose: print "OneD moments: h: %g  a: %g  c: %g  w: %g" % tuple(params)

        parinfo = [ {'n':0,'value':params[0],'limits':[minpars[0],maxpars[0]],'limited':[limitedmin[0],limitedmax[0]],'fixed':fixed[0],'parname':"HEIGHT",'error':0} ,
                    {'n':1,'value':params[1],'limits':[minpars[1],maxpars[1]],'limited':[limitedmin[1],limitedmax[1]],'fixed':fixed[1],'parname':"AMPLITUDE",'error':0},
                    {'n':2,'value':params[2],'limits':[minpars[2],maxpars[2]],'limited':[limitedmin[2],limitedmax[2]],'fixed':fixed[2],'parname':"SHIFT",'error':0},
                    {'n':3,'value':params[3],'limits':[minpars[3],maxpars[3]],'limited':[limitedmin[3],limitedmax[3]],'fixed':fixed[3],'parname':"GWIDTH",'error':0},
                    {'n':4,'value':params[4],'limits':[minpars[4],maxpars[4]],'limited':[limitedmin[4],limitedmax[4]],'fixed':fixed[4],'parname':"LWIDTH",'error':0}]

        mp = mpfit(mpfitfun(xax,data,err),parinfo=parinfo,quiet=quiet)
        mpp = mp.params
        if mp.perror is not None: mpperr = mp.perror
        else: mpperr = mpp*0
        chi2 = mp.fnorm

        if mp.status == 0:
            raise Exception(mp.errmsg)

        if (not shh) or veryverbose:
            print "Fit message ",mp.errmsg
            print "Fit status: ",mp.status
            for i,p in enumerate(mpp):
                parinfo[i]['value'] = p
                print parinfo[i]['parname'],p," +/- ",mpperr[i]
            print "Chi2: ",mp.fnorm," Reduced Chi2: ",mp.fnorm/len(data)," DOF:",len(data)-len(mpp)

        self.mp = mp
        self.mpp = mpp[1:]
        self.mpperr = mpperr
        self.model = self.n_voigt(pars=mpp[1:])(xax) 
        return mpp,self.n_voigt(pars=mpp[1:])(xax),mpperr,chi2


    def multivoigtfit(self,xax, data, npeaks=1, err=None, params=[1,0,1,1],
            fixed=[False,False,False,False],
            limitedmin=[False,False,True,True],
            limitedmax=[False,False,False,False], minpars=[0,0,0,0],
            maxpars=[0,0,0,0], quiet=True, shh=True, veryverbose=False,
            **kwargs):
        """
        Fit multiple voigt profiles

        Inputs:
           xax - x axis
           data - y axis
           npeaks - How many voigt profiles to fit?  Default 1 (this could supersede onedgaussfit)
           err - error corresponding to data

         These parameters need to have length = 3*npeaks.  If npeaks > 1 and length = 3, they will
         be replicated npeaks times, otherwise they will be reset to defaults:
           params - Fit parameters: [amplitude, offset, Gfwhm, Lfwhm] * npeaks
                  If len(params) % 4 == 0, npeaks will be set to len(params) / 3
           fixed - Is parameter fixed?
           limitedmin/minpars - set lower limits on each parameter (default: width>0)
           limitedmax/maxpars - set upper limits on each parameter

           quiet - should MPFIT output each iteration?
           shh - output final parameters?

        Returns:
           Fit parameters
           Model
           Fit errors
           chi2
        """

        if len(params) != npeaks and (len(params) / 4) > npeaks:
            npeaks = len(params) / 4 
        self.npeaks = npeaks

        if isinstance(params,numpy.ndarray): params=params.tolist()

        # make sure all various things are the right length; if they're not, fix them using the defaults
        for parlist in (params,fixed,limitedmin,limitedmax,minpars,maxpars):
            if len(parlist) != 4*npeaks:
                # if you leave the defaults, or enter something that can be multiplied by 3 to get to the
                # right number of gaussians, it will just replicate
                if len(parlist) == 4: 
                    parlist *= npeaks 
                elif parlist==params:
                    parlist[:] = [1,0,1,1] * npeaks
                elif parlist==fixed or parlist==limitedmax:
                    parlist[:] = [False,False,False,False] * npeaks
                elif parlist==limitedmin:
                    parlist[:] = [False,False,True,True] * npeaks
                elif parlist==minpars or parlist==maxpars:
                    parlist[:] = [0,0,0,0] * npeaks

        def mpfitfun(x,y,err):
            if err is None:
                def f(p,fjac=None): return [0,(y-self.n_voigt(pars=p)(x))]
            else:
                def f(p,fjac=None): return [0,(y-self.n_voigt(pars=p)(x))/err]
            return f

        if xax == None:
            xax = numpy.arange(len(data))

        parnames = {0:"AMPLITUDE",1:"SHIFT",2:"GFWHM",3:"LFWHM"}

        parinfo = [ {'n':ii, 'value':params[ii],
            'limits':[minpars[ii],maxpars[ii]],
            'limited':[limitedmin[ii],limitedmax[ii]], 'fixed':fixed[ii],
            'parname':parnames[ii%4]+str(ii/4), 'error':ii} 
            for ii in xrange(len(params)) ]

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
            print "Fit message ",mp.errmsg
            print "Final fit values: "
            for i,p in enumerate(mpp):
                parinfo[i]['value'] = p
                print parinfo[i]['parname'],p," +/- ",mpperr[i]
            print "Chi2: ",mp.fnorm," Reduced Chi2: ",mp.fnorm/len(data)," DOF:",len(data)-len(mpp)

        self.mp = mp
        self.mpp = mpp
        self.mpperr = mpperr
        self.model = self.n_voigt(pars=mpp)(xax) 
        return mpp,self.n_voigt(pars=mpp)(xax),mpperr,chi2

    def annotations(self):
        label_list = [(
                "$A(%i)$=%6.4g $\\pm$ %6.4g" % (jj,self.mpp[0+jj*self.npars],self.mpperr[0+jj*self.npars]),
                "$x(%i)$=%6.4g $\\pm$ %6.4g" % (jj,self.mpp[1+jj*self.npars],self.mpperr[1+jj*self.npars]),
                "$\\sigma_G(%i)$=%6.4g $\\pm$ %6.4g" % (jj,self.mpp[2+jj*self.npars],self.mpperr[2+jj*self.npars]),
                "$\\sigma_L(%i)$=%6.4g $\\pm$ %6.4g" % (jj,self.mpp[3+jj*self.npars],self.mpperr[3+jj*self.npars]),
                          ) for jj in range(self.npeaks)]
        labels = tuple(mpcb.flatten(label_list))
        return labels

    def components(self,xarr,modelpars):

        modelcomponents = [ self.voigt(xarr,
            modelpars[4*i],modelpars[4*i+1],modelpars[4*i+2],modelpars[4*i+3])
            for i in range(self.npeaks)]

        return modelcomponents
