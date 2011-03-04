# gaussfitter.py
# created by Adam Ginsburg (adam.ginsburg@colorado.edu or keflavich@gmail.com) 3/17/08)
# latest version available at http://code.google.com/p/agpy/source/browse/trunk/agpy/gaussfitter.py
import numpy
from numpy.ma import median
from numpy import pi
#from scipy import optimize,stats,pi
from mpfit import mpfit
""" 
Note about mpfit/leastsq: 
I switched everything over to the Markwardt mpfit routine for a few reasons,
but foremost being the ability to set limits on parameters, not just force them
to be fixed.  As far as I can tell, leastsq does not have that capability.  

The version of mpfit I use can be found here:
    http://code.google.com/p/agpy/source/browse/trunk/mpfit
"""

"""
To do:
    -turn into a class instead of a collection of objects
    -implement WCS-based gaussian fitting with correct coordinates
"""

def moments(data,circle,rotate,vheight,estimator=median,**kwargs):
    """Returns (height, amplitude, x, y, width_x, width_y, rotation angle)
    the gaussian parameters of a 2D distribution by calculating its
    moments.  Depending on the input parameters, will only output 
    a subset of the above.
    
    If using masked arrays, pass estimator=numpy.ma.median
    """
    total = numpy.abs(data).sum()
    Y, X = numpy.indices(data.shape) # python convention: reverse x,y numpy.indices
    y = numpy.argmax((X*numpy.abs(data)).sum(axis=1)/total)
    x = numpy.argmax((Y*numpy.abs(data)).sum(axis=0)/total)
    col = data[int(y),:]
    # FIRST moment, not second!
    width_x = numpy.sqrt(numpy.abs((numpy.arange(col.size)-y)*col).sum()/numpy.abs(col).sum())
    row = data[:, int(x)]
    width_y = numpy.sqrt(numpy.abs((numpy.arange(row.size)-x)*row).sum()/numpy.abs(row).sum())
    width = ( width_x + width_y ) / 2.
    height = estimator(data.ravel())
    amplitude = data.max()-height
    mylist = [amplitude,x,y]
    if numpy.isnan(width_y) or numpy.isnan(width_x) or numpy.isnan(height) or numpy.isnan(amplitude):
        raise ValueError("something is nan")
    if vheight==1:
        mylist = [height] + mylist
    if circle==0:
        mylist = mylist + [width_x,width_y]
        if rotate==1:
            mylist = mylist + [0.] #rotation "moment" is just zero...
            # also, circles don't rotate.
    else:  
        mylist = mylist + [width]
    return mylist

def twodgaussian(inpars, circle=0, rotate=1, vheight=1, shape=None):
    """Returns a 2d gaussian function of the form:
        x' = numpy.cos(rota) * x - numpy.sin(rota) * y
        y' = numpy.sin(rota) * x + numpy.cos(rota) * y
        (rota should be in degrees)
        g = b + a * numpy.exp ( - ( ((x-center_x)/width_x)**2 +
        ((y-center_y)/width_y)**2 ) / 2 )

        inpars = [b,a,center_x,center_y,width_x,width_y,rota]
                 (b is background height, a is peak amplitude)

        where x and y are the input parameters of the returned function,
        and all other parameters are specified by this function

        However, the above values are passed by list.  The list should be:
        inpars = (height,amplitude,center_x,center_y,width_x,width_y,rota)

        You can choose to ignore / neglect some of the above input parameters 
            unumpy.sing the following options:
            circle=0 - default is an elliptical gaussian (different x, y
                widths), but can reduce the input by one parameter if it's a
                circular gaussian
            rotate=1 - default allows rotation of the gaussian ellipse.  Can
                remove last parameter by setting rotate=0
            vheight=1 - default allows a variable height-above-zero, i.e. an
                additive constant for the Gaussian function.  Can remove first
                parameter by setting this to 0
            shape=None - if shape is set (to a 2-parameter list) then returns
                an image with the gaussian defined by inpars
        """
    inpars_old = inpars
    inpars = list(inpars)
    if vheight == 1:
        height = inpars.pop(0)
        height = float(height)
    else:
        height = float(0)
    amplitude, center_y, center_x = inpars.pop(0),inpars.pop(0),inpars.pop(0)
    amplitude = float(amplitude)
    center_x = float(center_x)
    center_y = float(center_y)
    if circle == 1:
        width = inpars.pop(0)
        width_x = float(width)
        width_y = float(width)
        rotate = 0
    else:
        width_x, width_y = inpars.pop(0),inpars.pop(0)
        width_x = float(width_x)
        width_y = float(width_y)
    if rotate == 1:
        rota = inpars.pop(0)
        rota = pi/180. * float(rota)
        rcen_x = center_x * numpy.cos(rota) - center_y * numpy.sin(rota)
        rcen_y = center_x * numpy.sin(rota) + center_y * numpy.cos(rota)
    else:
        rcen_x = center_x
        rcen_y = center_y
    if len(inpars) > 0:
        raise ValueError("There are still input parameters:" + str(inpars) + \
                " and you've input: " + str(inpars_old) + \
                " circle=%d, rotate=%d, vheight=%d" % (circle,rotate,vheight) )
            
    def rotgauss(x,y):
        if rotate==1:
            xp = x * numpy.cos(rota) - y * numpy.sin(rota)
            yp = x * numpy.sin(rota) + y * numpy.cos(rota)
        else:
            xp = x
            yp = y
        g = height+amplitude*numpy.exp(
            -(((rcen_x-xp)/width_x)**2+
            ((rcen_y-yp)/width_y)**2)/2.)
        return g
    if shape is not None:
        return rotgauss(*numpy.indices(shape))
    else:
        return rotgauss

def gaussfit(data,err=None,params=[],autoderiv=1,return_all=0,circle=0,
        fixed=numpy.repeat(False,7),limitedmin=[False,False,False,False,True,True,True],
        limitedmax=[False,False,False,False,False,False,True],
        usemoment=numpy.array([],dtype='bool'),
        minpars=numpy.repeat(0,7),maxpars=[0,0,0,0,0,0,360],
        rotate=1,vheight=1,quiet=True,returnmp=False,
        returnfitimage=False,**kwargs):
    """
    Gaussian fitter with the ability to fit a variety of different forms of
    2-dimensional gaussian.
    
    Input Parameters:
        data - 2-dimensional data array
        err=None - error array with same size as data array
        params=[] - initial input parameters for Gaussian function.
            (height, amplitude, x, y, width_x, width_y, rota)
            if not input, these will be determined from the moments of the system, 
            assuming no rotation
        autoderiv=1 - use the autoderiv provided in the lmder.f function (the
            alternative is to us an analytic derivative with lmdif.f: this method
            is less robust)
        return_all=0 - Default is to return only the Gaussian parameters.  
                   1 - fit params, fit error
        returnfitimage - returns (best fit params,best fit image)
        returnmp - returns the full mpfit struct
        circle=0 - default is an elliptical gaussian (different x, y widths),
            but can reduce the input by one parameter if it's a circular gaussian
        rotate=1 - default allows rotation of the gaussian ellipse.  Can remove
            last parameter by setting rotate=0.  numpy.expects angle in DEGREES
        vheight=1 - default allows a variable height-above-zero, i.e. an
            additive constant for the Gaussian function.  Can remove first
            parameter by setting this to 0
        usemoment - can choose which parameters to use a moment estimation for.
            Other parameters will be taken from params.  Needs to be a boolean
            array.

    Output:
        Default output is a set of Gaussian parameters with the same shape as
            the input parameters

        Can also output the covariance matrix, 'infodict' that contains a lot
            more detail about the fit (see scipy.optimize.leastsq), and a message
            from leastsq telling what the exit status of the fitting routine was

        Warning: Does NOT necessarily output a rotation angle between 0 and 360 degrees.
    """
    usemoment=numpy.array(usemoment,dtype='bool')
    params=numpy.array(params,dtype='float')
    if usemoment.any() and len(params)==len(usemoment):
        moment = numpy.array(moments(data,circle,rotate,vheight,**kwargs),dtype='float')
        params[usemoment] = moment[usemoment]
    elif params == [] or len(params)==0:
        params = (moments(data,circle,rotate,vheight,**kwargs))
    if vheight==0:
        vheight=1
        params = numpy.concatenate([[0],params])
        fixed[0] = 1


    # mpfit will fail if it is given a start parameter outside the allowed range:
    for i in xrange(len(params)): 
        if params[i] > maxpars[i] and limitedmax[i]: params[i] = maxpars[i]
        if params[i] < minpars[i] and limitedmin[i]: params[i] = minpars[i]

    if err == None:
        errorfunction = lambda p: numpy.ravel((twodgaussian(p,circle,rotate,vheight)\
                (*numpy.indices(data.shape)) - data))
    else:
        errorfunction = lambda p: numpy.ravel((twodgaussian(p,circle,rotate,vheight)\
                (*numpy.indices(data.shape)) - data)/err)
    def mpfitfun(data,err):
        if err == None:
            def f(p,fjac=None): return [0,numpy.ravel(data-twodgaussian(p,circle,rotate,vheight)\
                    (*numpy.indices(data.shape)))]
        else:
            def f(p,fjac=None): return [0,numpy.ravel((data-twodgaussian(p,circle,rotate,vheight)\
                    (*numpy.indices(data.shape)))/err)]
        return f

                    
    parinfo = [ 
                {'n':1,'value':params[1],'limits':[minpars[1],maxpars[1]],'limited':[limitedmin[1],limitedmax[1]],'fixed':fixed[1],'parname':"AMPLITUDE",'error':0},
                {'n':2,'value':params[2],'limits':[minpars[2],maxpars[2]],'limited':[limitedmin[2],limitedmax[2]],'fixed':fixed[2],'parname':"XSHIFT",'error':0},
                {'n':3,'value':params[3],'limits':[minpars[3],maxpars[3]],'limited':[limitedmin[3],limitedmax[3]],'fixed':fixed[3],'parname':"YSHIFT",'error':0},
                {'n':4,'value':params[4],'limits':[minpars[4],maxpars[4]],'limited':[limitedmin[4],limitedmax[4]],'fixed':fixed[4],'parname':"XWIDTH",'error':0} ]
    if vheight == 1:
        parinfo.insert(0,{'n':0,'value':params[0],'limits':[minpars[0],maxpars[0]],'limited':[limitedmin[0],limitedmax[0]],'fixed':fixed[0],'parname':"HEIGHT",'error':0})
    if circle == 0:
        parinfo.append({'n':5,'value':params[5],'limits':[minpars[5],maxpars[5]],'limited':[limitedmin[5],limitedmax[5]],'fixed':fixed[5],'parname':"YWIDTH",'error':0})
        if rotate == 1:
            parinfo.append({'n':6,'value':params[6],'limits':[minpars[6],maxpars[6]],'limited':[limitedmin[6],limitedmax[6]],'fixed':fixed[6],'parname':"ROTATION",'error':0})

    if autoderiv == 0:
        # the analytic derivative, while not terribly difficult, is less
        # efficient and useful.  I only bothered putting it here because I was
        # instructed to do so for a class project - please ask if you would
        # like this feature implemented
        raise ValueError("I'm sorry, I haven't implemented this feature yet.")
    else:
#        p, cov, infodict, errmsg, success = optimize.leastsq(errorfunction,\
#                params, full_output=1)
        mp = mpfit(mpfitfun(data,err),parinfo=parinfo,quiet=quiet)


    if returnmp:
        returns = (mp)
    elif return_all == 0:
        returns = mp.params
    elif return_all == 1:
        returns = mp.params,mp.perror
    if returnfitimage:
        fitimage = twodgaussian(mp.params,circle,rotate,vheight)(*numpy.indices(data.shape))
        returns = (returns,fitimage)
    return returns

def onedmoments(Xax,data,vheight=True,estimator=median,negamp=None,
        veryverbose=False, **kwargs):
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

def onedgaussian(x,H,A,dx,w):
    """
    Returns a 1-dimensional gaussian of form
    H+A*numpy.exp(-(x-dx)**2/(2*w**2))
    """
    return H+A*numpy.exp(-(x-dx)**2/(2*w**2))

def onedgaussfit(xax, data, err=None,
        params=[0,1,0,1],fixed=[False,False,False,False],
        limitedmin=[False,False,False,True],
        limitedmax=[False,False,False,False], minpars=[0,0,0,0],
        maxpars=[0,0,0,0], quiet=True, shh=True,
        veryverbose=False,
        vheight=True, negamp=False,
        usemoments=False):
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
        if err == None:
            def f(p,fjac=None): return [0,(y-onedgaussian(x,*p))]
        else:
            def f(p,fjac=None): return [0,(y-onedgaussian(x,*p))/err]
        return f

    if xax == None:
        xax = numpy.arange(len(data))

    if vheight is False: 
        height = params[0]
        fixed[0] = True
    if usemoments:
        params = onedmoments(xax,data,vheight=vheight,negamp=negamp, veryverbose=veryverbose)
        if vheight is False: params = [height]+params
        if veryverbose: print "OneD moments: h: %g  a: %g  c: %g  w: %g" % tuple(params)

    parinfo = [ {'n':0,'value':params[0],'limits':[minpars[0],maxpars[0]],'limited':[limitedmin[0],limitedmax[0]],'fixed':fixed[0],'parname':"HEIGHT",'error':0} ,
                {'n':1,'value':params[1],'limits':[minpars[1],maxpars[1]],'limited':[limitedmin[1],limitedmax[1]],'fixed':fixed[1],'parname':"AMPLITUDE",'error':0},
                {'n':2,'value':params[2],'limits':[minpars[2],maxpars[2]],'limited':[limitedmin[2],limitedmax[2]],'fixed':fixed[2],'parname':"SHIFT",'error':0},
                {'n':3,'value':params[3],'limits':[minpars[3],maxpars[3]],'limited':[limitedmin[3],limitedmax[3]],'fixed':fixed[3],'parname':"WIDTH",'error':0}]

    mp = mpfit(mpfitfun(xax,data,err),parinfo=parinfo,quiet=quiet)
    mpp = mp.params
    mpperr = mp.perror
    chi2 = mp.fnorm

    if mp.status == 0:
        raise Exception(mp.errmsg)

    if (not shh) or veryverbose:
        print "Fit status: ",mp.status
        for i,p in enumerate(mpp):
            parinfo[i]['value'] = p
            print parinfo[i]['parname'],p," +/- ",mpperr[i]
        print "Chi2: ",mp.fnorm," Reduced Chi2: ",mp.fnorm/len(data)," DOF:",len(data)-len(mpp)

    return mpp,onedgaussian(xax,*mpp),mpperr,chi2


def n_gaussian(pars=None,a=None,dx=None,sigma=None):
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

def multigaussfit(xax, data, ngauss=1, err=None, params=[1,0,1],
        fixed=[False,False,False], limitedmin=[False,False,True],
        limitedmax=[False,False,False], minpars=[0,0,0], maxpars=[0,0,0],
        quiet=True, shh=True, veryverbose=False):
    """
    An improvement on onedgaussfit.  Lets you fit multiple gaussians.

    Inputs:
       xax - x axis
       data - y axis
       ngauss - How many gaussians to fit?  Default 1 (this could supersede onedgaussfit)
       err - error corresponding to data

     These parameters need to have length = 3*ngauss.  If ngauss > 1 and length = 3, they will
     be replicated ngauss times, otherwise they will be reset to defaults:
       params - Fit parameters: [amplitude, offset, width] * ngauss
              If len(params) % 3 == 0, ngauss will be set to len(params) / 3
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

    if len(params) != ngauss and (len(params) / 3) > ngauss:
        ngauss = len(params) / 3 

    if isinstance(params,numpy.ndarray): params=params.tolist()

    # make sure all various things are the right length; if they're not, fix them using the defaults
    for parlist in (params,fixed,limitedmin,limitedmax,minpars,maxpars):
        if len(parlist) != 3*ngauss:
            # if you leave the defaults, or enter something that can be multiplied by 3 to get to the
            # right number of gaussians, it will just replicate
            if len(parlist) == 3: 
                parlist *= ngauss 
            elif parlist==params:
                parlist[:] = [1,0,1] * ngauss
            elif parlist==fixed or parlist==limitedmax:
                parlist[:] = [False,False,False] * ngauss
            elif parlist==limitedmin:
                parlist[:] = [False,False,True] * ngauss
            elif parlist==minpars or parlist==maxpars:
                parlist[:] = [0,0,0] * ngauss

    def mpfitfun(x,y,err):
        if err == None:
            def f(p,fjac=None): return [0,(y-n_gaussian(pars=p)(x))]
        else:
            def f(p,fjac=None): return [0,(y-n_gaussian(pars=p)(x))/err]
        return f

    if xax == None:
        xax = numpy.arange(len(data))

    parnames = {0:"AMPLITUDE",1:"SHIFT",2:"WIDTH"}

    parinfo = [ {'n':ii, 'value':params[ii],
        'limits':[minpars[ii],maxpars[ii]],
        'limited':[limitedmin[ii],limitedmax[ii]], 'fixed':fixed[ii],
        'parname':parnames[ii%3]+str(ii%3), 'error':ii} 
        for ii in xrange(len(params)) ]

    if veryverbose:
        print "GUESSES: "
        print "\n".join(["%s: %s" % (p['parname'],p['value']) for p in parinfo])

    mp = mpfit(mpfitfun(xax,data,err),parinfo=parinfo,quiet=quiet)
    mpp = mp.params
    mpperr = mp.perror
    chi2 = mp.fnorm

    if mp.status == 0:
        raise Exception(mp.errmsg)

    if not shh:
        print "Final fit values: "
        for i,p in enumerate(mpp):
            parinfo[i]['value'] = p
            print parinfo[i]['parname'],p," +/- ",mpperr[i]
        print "Chi2: ",mp.fnorm," Reduced Chi2: ",mp.fnorm/len(data)," DOF:",len(data)-len(mpp)

    return mpp,n_gaussian(pars=mpp)(xax),mpperr,chi2

