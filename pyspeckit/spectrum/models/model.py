import numpy as np

class SpectralModel(object):

    def __init__(self, modelfunc, npars, parnames=None, parvalues=None, parlimits=None,
            parlimited=None, parfixed=None, parerror=None):

        self.modelfunc = modelfunc
        self.npars = npars
        self.parnames = parnames

        for var in (parnames,parvalues,parlimits,parlimited,parfixed,parerror):
            if var is None:
                var = np.zeros(self.npars, dtype='bool')

        self.parinfo = [ {'n':ii,
            'value':parvalues[ii],
            'limits':parlimits[ii],
            'limited':parlimited[ii],
            'fixed':parfixed[ii],
            'parname':parnames[ii],
            'error':parerror[ii]} 
            for ii in xrange(npars)]

    def mpfitfun(x,y,err):
        if err is None:
            def f(p,fjac=None): return [0,(y-self.modelfunc(x,*p))]
        else:
            def f(p,fjac=None): return [0,(y-self.modelfunc(x,*p))/err]
        return f

    def __call__(self, xax, data, quiet=True, shh=True, veryverbose=False, **kwargs):

        mp = mpfit(self.mpfitfun(xax,data,err),parinfo=self.parinfo,quiet=quiet,**kwargs)
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
                self.parinfo[i]['value'] = p
                print self.parinfo[i]['parname'],p," +/- ",mpperr[i]
            print "Chi2: ",mp.fnorm," Reduced Chi2: ",mp.fnorm/len(data)," DOF:",len(data)-len(mpp)

        self.mp = mp
        self.mpp = mpp[1:]
        self.mpperr = mpperr[1:]
        self.model = self.modelfunc(xax,*mpp)
        return mpp,self.model,mpperr,chi2


