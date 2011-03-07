import numpy
from numpy.ma import median
from numpy import pi
from mpfit import mpfit

try:
    import scipy.special
except ImportError:
    print "Voigt profile fitting will fail.  Requires scipy.special"


def voigt(xarr,amp,xcen,Gfwhm,Lfwhm):
    """
    voigt profile

    V(x,sig,gam) = Re(w(z))/(sig*sqrt(2*pi))
    z = (x+i*gam)/(sig*sqrt(2))

    Converted from 
    http://mail.scipy.org/pipermail/scipy-user/2011-January/028327.html
    """

    tmp = 1.0/scipy.special.wofz(numpy.zeros((len(xarr))) \
          +1j*numpy.sqrt(numpy.log(2.0))*Lfwhm).real
    tmp = tmp*amp* \
          scipy.special.wofz(2*numpy.sqrt(numpy.log(2.0))*(xarr-xcen)/Gfwhm+1j* \
          numpy.sqrt(numpy.log(2.0))*Lfwhm).real
    return tmp

def n_voigt(pars=None):
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
            v += voigt(x,a[i],dx[i],Gfwhm[i],Lfwhm[i])
        return v
    return L

def multivoigtfit(xax, data, nvoigt=1, err=None, params=[1,0,1,1],
        fixed=[False,False,False,False],
        limitedmin=[False,False,True,True],
        limitedmax=[False,False,False,False], minpars=[0,0,0,0],
        maxpars=[0,0,0,0], quiet=True, shh=True, veryverbose=False):
    """
    Fit multiple voigt profiles

    Inputs:
       xax - x axis
       data - y axis
       nvoigt - How many voigt profiles to fit?  Default 1 (this could supersede onedgaussfit)
       err - error corresponding to data

     These parameters need to have length = 3*nvoigt.  If nvoigt > 1 and length = 3, they will
     be replicated nvoigt times, otherwise they will be reset to defaults:
       params - Fit parameters: [amplitude, offset, Gfwhm, Lfwhm] * nvoigt
              If len(params) % 4 == 0, nvoigt will be set to len(params) / 3
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

    if len(params) != nvoigt and (len(params) / 4) > nvoigt:
        nvoigt = len(params) / 4 

    if isinstance(params,numpy.ndarray): params=params.tolist()

    # make sure all various things are the right length; if they're not, fix them using the defaults
    for parlist in (params,fixed,limitedmin,limitedmax,minpars,maxpars):
        if len(parlist) != 4*nvoigt:
            # if you leave the defaults, or enter something that can be multiplied by 3 to get to the
            # right number of gaussians, it will just replicate
            if len(parlist) == 4: 
                parlist *= nvoigt 
            elif parlist==params:
                parlist[:] = [1,0,1,1] * nvoigt
            elif parlist==fixed or parlist==limitedmax:
                parlist[:] = [False,False,False,False] * nvoigt
            elif parlist==limitedmin:
                parlist[:] = [False,False,True,True] * nvoigt
            elif parlist==minpars or parlist==maxpars:
                parlist[:] = [0,0,0,0] * nvoigt

    def mpfitfun(x,y,err):
        if err == None:
            def f(p,fjac=None): return [0,(y-n_voigt(pars=p)(x))]
        else:
            def f(p,fjac=None): return [0,(y-n_voigt(pars=p)(x))/err]
        return f

    if xax == None:
        xax = numpy.arange(len(data))

    parnames = {0:"AMPLITUDE",1:"SHIFT",2:"GFWHM",3:"LFWHM"}

    parinfo = [ {'n':ii, 'value':params[ii],
        'limits':[minpars[ii],maxpars[ii]],
        'limited':[limitedmin[ii],limitedmax[ii]], 'fixed':fixed[ii],
        'parname':parnames[ii%4]+str(ii%4), 'error':ii} 
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
        print "Final fit values: "
        for i,p in enumerate(mpp):
            parinfo[i]['value'] = p
            print parinfo[i]['parname'],p," +/- ",mpperr[i]
        print "Chi2: ",mp.fnorm," Reduced Chi2: ",mp.fnorm/len(data)," DOF:",len(data)-len(mpp)

    return mpp,n_voigt(pars=mpp)(xax),mpperr,chi2

