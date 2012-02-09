"""
===============
Gaussian fitter
===============

.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>
Created 3/17/08

Original version available at http://code.google.com/p/agpy/source/browse/trunk/agpy/gaussfitter.py
(the version below uses a Class instead of independent functions)
"""
import numpy
from numpy.ma import median
from numpy import pi
from pyspeckit.mpfit import mpfit
import matplotlib.cbook as mpcb
from . import mpfit_messages
from . import model

class gaussian_fitter(model.SpectralModel):
    """
    A rather complicated Gaussian fitter class.  Inherits from, but overrides
    most components of, :mod:`model.SpectralModel`
    """

    def __init__(self,multisingle='multi'):
        self.npars = 3
        self.npeaks = 1
        self.onepeakgaussfit = self._fourparfitter(self.onepeakgaussian)
        if multisingle in ('multi','single'):
            self.multisingle = multisingle
        else:
            raise Exception("multisingle must be multi or single")

    def __call__(self,*args,**kwargs):
        if self.multisingle == 'single':
            return self.onepeakgaussfit(*args,**kwargs)
        elif self.multisingle == 'multi':
            return self.multigaussfit(*args,**kwargs)

    def onepeakgaussian(self, x,H,A,dx,w):
        """
        Returns a 1-dimensional gaussian of form
        H+A*numpy.exp(-(x-dx)**2/(2*w**2))
        
        [height,amplitude,center,width]
        
        """
        x = numpy.array(x) # make sure xarr is no longer a spectroscopic axis
        return H+A*numpy.exp(-(x-dx)**2/(2*w**2))
        
    def multipeakgaussian(self, x, pars):
        """
        Returns flux at position x due to contributions from multiple Gaussians.
        """
        x = numpy.array(x) # make sure xarr is no longer a spectroscopic axis
        
        pars = numpy.reshape(pars, (len(pars) / 3, 3))
        
        result = 0
        for fit in pars: result += self.onepeakgaussian(x, 0, fit[0], fit[1], fit[2])
        return result
        
    def slope(self, x):
        """
        Return slope at position x for multicomponent Gaussian fit.  Need this in measurements class for
        finding the FWHM of multicomponent lines whose centroids are not identical.
        """    
        
        pars = numpy.reshape(self.mpp, (len(self.mpp) / 3, 3))
        result = 0
        for fit in pars:
            result += self.onepeakgaussian(x, 0, fit[0], fit[1], fit[2]) * (-2. * (x - fit[1]) / 2. / fit[2]**2)
        return result

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
            for ii in range(len(pars)/3):
                v += a[ii] * numpy.exp( - ( x - dx[ii] )**2 / (2.0*sigma[ii]**2) )
            return v
        return g

    def multigaussfit(self, xax, data, npeaks=1, err=None, params=[1,0,1],
            fixed=[False,False,False], limitedmin=[False,False,True],
            limitedmax=[False,False,False], minpars=[0,0,0], maxpars=[0,0,0],
            quiet=True, shh=True, veryverbose=False, negamp=None,
            tied = ['', '', ''], parinfo=None, debug=False, **kwargs):
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
        
        if len(params) != npeaks and (len(params) / 3) > npeaks:
            self.npeaks = len(params) / 3 
        else:
            self.npeaks = npeaks

        if isinstance(params,numpy.ndarray): params=params.tolist()

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
                def f(p,fjac=None): return [0,(y-self.n_gaussian(pars=p)(x))]
            else:
                def f(p,fjac=None): return [0,(y-self.n_gaussian(pars=p)(x))/err]
            return f

        if xax is None:
            xax = numpy.arange(len(data))

        parnames = {0:"AMPLITUDE",1:"SHIFT",2:"WIDTH"}

        if parinfo is None:
            parinfo = [ {'n':ii, 'value':params[ii],
                'limits':[minpars[ii],maxpars[ii]],
                'limited':[limitedmin[ii],limitedmax[ii]], 'fixed':fixed[ii],
                'parname':parnames[ii%3]+str(ii/3), 'error':ii, 'tied':tied[ii]} 
                for ii in xrange(len(params)) ]

        if veryverbose:
            print "GUESSES: "
            print "\n".join(["%s: %s" % (p['parname'],p['value']) for p in parinfo])

        if debug: 
            for p in parinfo: print p
            
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
                "$x(%i)$=%6.4g $\\pm$ %6.4g" % (jj,self.mpp[1+jj*self.npars],self.mpperr[1+jj*self.npars]),
                "$\\sigma(%i)$=%6.4g $\\pm$ %6.4g" % (jj,self.mpp[2+jj*self.npars],self.mpperr[2+jj*self.npars])
                          ) for jj in range(self.npeaks)]
        labels = tuple(mpcb.flatten(label_list))
        return labels

    def components(self,xarr,modelpars):

        modelcomponents = [ self.onepeakgaussian(xarr,
            0.0,modelpars[3*i],modelpars[3*i+1],modelpars[3*i+2]) for i in range(self.npeaks)]

        return modelcomponents

    def integral(self, modelpars):
        """
        Return the integral of the individual components (ignoring height)
        """
        return self.model.sum()

        # this is the "proper" way to do it, but the above line was used for compatibility with other models
        integ = 0
        if len(modelpars) % 3 == 0:
            for amp,cen,width in numpy.reshape(modelpars,[len(modelpars)/3,3]):
                integ += amp*width*numpy.sqrt(2.0*numpy.pi)

        return integ

    n_modelfunc = n_gaussian

