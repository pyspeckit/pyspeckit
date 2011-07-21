import numpy
from mpfit import mpfit
from numpy.ma import median

class SimpleFitter(object):

    def __init__():
        pass

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
        if negamp and amplitude > 0 and veryverbose: 
            print "WARNING: likely fit failure.  negamp=True, but amplitude > 0"
        if negamp is False and amplitude < 0 and veryverbose: 
            print "WARNING: likely fit failure.  negamp=False, but amplitude < 0"
        if numpy.isnan(width_x) or numpy.isnan(height) or numpy.isnan(amplitude):
            raise ValueError("something is nan")
        if vheight:
            mylist = [height] + mylist
        return mylist

    def _fourparfitter(self, modelfunc):

        def fitter(xax, data, err=None,
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
                    def f(p,fjac=None): return [0,(y-modelfunc(x,*p))]
                else:
                    def f(p,fjac=None): return [0,(y-modelfunc(x,*p))/err]
                return f

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
            self.model = modelfunc(xax,*mpp)
            return mpp,modelfunc(xax,*mpp),mpperr,chi2

        return fitter

def vheightmodel(zeroheightmodel):
    def vhm(xax, *pars,**kwargs):
        """
        Wrapper function vhm to set variable height.
        Parameter order: height, amplitude, shift, width
        """
        vheight=True
        if 'vheight' in kwargs:
            vheight = kwargs.pop('vheight')
        if vheight:
            return zeroheightmodel(xax, *pars[1:],**kwargs) + pars[0]
        else:
            return zeroheightmodel(xax, *pars[1:],**kwargs)
    vhm.__doc__ += zeroheightmodel.__doc__
    return vhm
