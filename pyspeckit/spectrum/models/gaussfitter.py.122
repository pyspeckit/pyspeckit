# gaussfitter.py
# created by Adam Ginsburg (adam.ginsburg@colorado.edu or keflavich@gmail.com) 3/17/08)
# latest version available at http://code.google.com/p/agpy/source/browse/trunk/agpy/gaussfitter.py
# (the version below uses a Class instead of independent functions)
import numpy
from numpy.ma import median
from numpy import pi
from mpfit import mpfit
import matplotlib.cbook as mpcb
from . import mpfit_messages

class gaussian_fitter(object):

    def __init__(self,multisingle='multi'):
        self.npars = 3
        self.npeaks = 1
        if multisingle in ('multi','single'):
            self.multisingle = multisingle
        else:
            raise Exception("multisingle must be multi or single")

    def __call__(self,*args,**kwargs):
        if self.multisingle == 'single':
            return self.onepeakgaussfit(*args,**kwargs)
        elif self.multisingle == 'multi':
            return self.multigaussfit(*args,**kwargs)

    def moments(self, Xax, data, vheight=True, estimator=median, negamp=None,
            veryverbose=False,  **kwargs):
        """Returns (height, amplitude, x, width_x)
        the gaussian parameters of a 1D distribution by calculating its
        moments.  Depending on the input parameters, will only output 
        a subset of the above.
        
        If using masked arrays, pass estimator=numpy.ma.median
        'estimator' is used to measure the background level (height)

        negamp can be used to force the peak negative (True), positive (False),
        or it will be "autodetected" (negamp=None)
        """

        dx = numpy.mean(Xax[1:] - Xax[:-1]) # assume a regular grid
        integral = (data*dx).sum()
        height = estimator(data)
        
        # try to figure out whether pos or neg based on the minimum width of the pos/neg peaks
        Lpeakintegral = integral - height*len(Xax)*dx - (data[data>height]*dx).sum()
        Lamplitude = data.min()-height
        Lwidth_x = 0.5*(numpy.abs(Lpeakintegral / Lamplitude))
        Hpeakintegral = integral - height*len(Xax)*dx - (data[data<height]*dx).sum()
        Hamplitude = data.max()-height
        Hwidth_x = 0.5*(numpy.abs(Hpeakintegral / Hamplitude))
        Lstddev = Xax[data<data.mean()].std()
        Hstddev = Xax[data>data.mean()].std()
        #print "Lstddev: %10.3g  Hstddev: %10.3g" % (Lstddev,Hstddev)
        #print "Lwidth_x: %10.3g  Hwidth_x: %10.3g" % (Lwidth_x,Hwidth_x)

        if negamp: # can force the guess to be negative
            xcen,amplitude,width_x = Xax[numpy.argmin(data)],Lamplitude,Lwidth_x
        elif negamp is None:
            if Hstddev < Lstddev: 
                xcen,amplitude,width_x, = Xax[numpy.argmax(data)],Hamplitude,Hwidth_x
            else:                                                                   
                xcen,amplitude,width_x, = Xax[numpy.argmin(data)],Lamplitude,Lwidth_x
        else:  # if negamp==False, make positive
            xcen,amplitude,width_x = Xax[numpy.argmax(data)],Hamplitude,Hwidth_x

        if veryverbose:
            print "negamp: %s  amp,width,cen Lower: %g, %g   Upper: %g, %g  Center: %g" %\
                    (negamp,Lamplitude,Lwidth_x,Hamplitude,Hwidth_x,xcen)
        mylist = [amplitude,xcen,width_x]
        if numpy.isnan(width_x) or numpy.isnan(height) or numpy.isnan(amplitude):
            raise ValueError("something is nan")
        if vheight:
            mylist = [height] + mylist
        return mylist

    def onepeakgaussian(self, x,H,A,dx,w):
        """
        Returns a 1-dimensional gaussian of form
        H+A*numpy.exp(-(x-dx)**2/(2*w**2))
        
        [height,amplitude,center,width]
        
        """
        return H+A*numpy.exp(-(x-dx)**2/(2*w**2))

    def onepeakgaussfit(self, xax, data, err=None,
            params=[0,1,0,1],fixed=[False,False,False,False],
            limitedmin=[False,False,False,True],
            limitedmax=[False,False,False,False], minpars=[0,0,0,0],
            maxpars=[0,0,0,0], quiet=True, shh=True,
            veryverbose=False,
            vheight=True, negamp=False,
            usemoments=False,**kwargs):
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

           kwargs are passed to mpfit

        Returns:
           Fit parameters
           Model
           Fit errors
           chi2
        """

        def mpfitfun(x,y,err):
            if err is None:
                def f(p,fjac=None): return [0,(y-self.onepeakgaussian(x,*p))]
            else:
                def f(p,fjac=None): return [0,(y-self.onepeakgaussian(x,*p))/err]
            return f

        if xax is None:
            xax = numpy.arange(len(data))

        if vheight is False: 
            height = params[0]
            fixed[0] = True
        if usemoments:
            params = moments(xax,data,vheight=vheight,negamp=negamp, veryverbose=veryverbose)
            if vheight is False: params = [height]+params
            if veryverbose: print "OneD moments: h: %g  a: %g  c: %g  w: %g" % tuple(params)

        parinfo = [ {'n':0,'value':params[0],'limits':[minpars[0],maxpars[0]],'limited':[limitedmin[0],limitedmax[0]],'fixed':fixed[0],'parname':"HEIGHT",'error':0} ,
                    {'n':1,'value':params[1],'limits':[minpars[1],maxpars[1]],'limited':[limitedmin[1],limitedmax[1]],'fixed':fixed[1],'parname':"AMPLITUDE",'error':0},
                    {'n':2,'value':params[2],'limits':[minpars[2],maxpars[2]],'limited':[limitedmin[2],limitedmax[2]],'fixed':fixed[2],'parname':"SHIFT",'error':0},
                    {'n':3,'value':params[3],'limits':[minpars[3],maxpars[3]],'limited':[limitedmin[3],limitedmax[3]],'fixed':fixed[3],'parname':"WIDTH",'error':0}]

        mp = mpfit(mpfitfun(xax,data,err),parinfo=parinfo,quiet=quiet,**kwargs)
        mpp = mp.params
        if mp.perror is not None: mpperr = mp.perror
        else: mpperr = mpp*0
        chi2 = mp.fnorm

        if mp.status == 0:
            raise Exception(mp.errmsg)

        if (not shh) or veryverbose:
            print "Fit status: ",mp.status
            print "Fit error message: ",mp.errmsg
            print "Fit message: ",mpfit_messages[mp.status]
            for i,p in enumerate(mpp):
                parinfo[i]['value'] = p
                print parinfo[i]['parname'],p," +/- ",mpperr[i]
            print "Chi2: ",mp.fnorm," Reduced Chi2: ",mp.fnorm/len(data)," DOF:",len(data)-len(mpp)

        self.mp = mp
        self.mpp = mpp[1:]
        self.mpperr = mpperr[1:]
        self.model = self.onepeakgaussian(xax,*mpp)
        return mpp,self.onepeakgaussian(xax,*mpp),mpperr,chi2


    def n_gaussian(self, pars=None,a=None,dx=None,sigma=None):
        """
        Returns a function that sums over N gaussians, where N is the length of
        a,dx,sigma *OR* N = len(pars) / 3

        The background "height" is assumed to be zero (you must "baseline" your
        spectrum before fitting)

        pars  - a list with len(pars) = 3n, assuming a,dx,sigma repeated
        dx    - offset (velocity center) values
        sigma - line widths
        a     - amplitudes
        """
        if len(pars) % 3 == 0:
            a = [pars[ii] for ii in xrange(0,len(pars),3)]
            dx = [pars[ii] for ii in xrange(1,len(pars),3)]
            sigma = [pars[ii] for ii in xrange(2,len(pars),3)]
        elif not(len(dx) == len(sigma) == len(a)):
            raise ValueError("Wrong array lengths! dx: %i  sigma: %i  a: %i" % (len(dx),len(sigma),len(a)))

        def g(x):
            v = numpy.zeros(len(x))
            for i in range(len(dx)):
                v += a[i] * numpy.exp( - ( x - dx[i] )**2 / (2.0*sigma[i]**2) )
            return v
        return g

    def multigaussfit(self, xax, data, npeaks=1, err=None, params=[1,0,1],
            fixed=[False,False,False], limitedmin=[False,False,True],
            limitedmax=[False,False,False], minpars=[0,0,0], maxpars=[0,0,0],
            quiet=True, shh=True, veryverbose=False, tied = ['', '', ''], **kwargs):
        """
        An improvement on onepeakgaussfit.  Lets you fit multiple gaussians.

        Inputs:
           xax - x axis
           data - y axis
           npeaks - How many gaussians to fit?  Default 1 (this could supersede onepeakgaussfit)
           err - error corresponding to data

         These parameters need to have length = 3*npeaks.  If npeaks > 1 and length = 3, they will
         be replicated npeaks times, otherwise they will be reset to defaults:
           params - Fit parameters: [amplitude, offset, width] * npeaks
                  If len(params) % 3 == 0, npeaks will be set to len(params) / 3
           fixed - Is parameter fixed?
           limitedmin/minpars - set lower limits on each parameter (default: width>0)
           limitedmax/maxpars - set upper limits on each parameter

           quiet - should MPFIT output each iteration?
           shh - output final parameters?

           kwargs are passed to mpfit

        Returns:
           Fit parameters
           Model
           Fit errors
           chi2
        """

        if len(params) != npeaks and (len(params) / 3) > npeaks:
            self.npeaks = len(params) / 3 
        else:
            self.npeaks = npeaks

        if isinstance(params,numpy.ndarray): params=params.tolist()

        # make sure all various things are the right length; if they're not, fix them using the defaults
        for parlist in (params,fixed,limitedmin,limitedmax,minpars,maxpars,tied):
            if len(parlist) != 3*self.npeaks:
                # if you leave the defaults, or enter something that can be multiplied by 3 to get to the
                # right number of gaussians, it will just replicate
                if len(parlist) == 3: 
                    parlist *= self.npeaks 
                # is any of this stuff valid?  I don't think so...
                elif parlist==params:
                    parlist[:] = [1,0,1] * self.npeaks
                elif parlist==fixed or parlist==limitedmax:
                    parlist[:] = [False,False,False] * self.npeaks
                elif parlist==limitedmin:
                    parlist[:] = [False,False,True] * self.npeaks
                elif parlist==minpars or parlist==maxpars:
                    parlist[:] = [0,0,0] * self.npeaks
                elif parlist==tied:
                    parlist[:] = ['','',''] * self.npeaks

        def mpfitfun(x,y,err):
            if err is None:
                def f(p,fjac=None): return [0,(y-self.n_gaussian(pars=p)(x))]
            else:
                def f(p,fjac=None): return [0,(y-self.n_gaussian(pars=p)(x))/err]
            return f

        if xax is None:
            xax = numpy.arange(len(data))

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
        self.model = self.n_gaussian(pars=mpp)(xax)
        return mpp,self.n_gaussian(pars=mpp)(xax),mpperr,chi2

    def annotations(self):
        label_list = [(
                "$A(%i)$=%6.4g $\\pm$ %6.4g" % (jj,self.mpp[0+jj*self.npars],self.mpperr[0+jj*self.npars]),
                "$v(%i)$=%6.4g $\\pm$ %6.4g" % (jj,self.mpp[1+jj*self.npars],self.mpperr[1+jj*self.npars]),
                "$\\sigma(%i)$=%6.4g $\\pm$ %6.4g" % (jj,self.mpp[2+jj*self.npars],self.mpperr[2+jj*self.npars])
                          ) for jj in range(self.npeaks)]
        labels = tuple(mpcb.flatten(label_list))
        return labels
