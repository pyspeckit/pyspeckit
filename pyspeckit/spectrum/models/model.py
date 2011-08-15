import numpy as np
from mpfit import mpfit
import copy

class SpectralModel(object):

    def __init__(self, modelfunc, npars, parnames=None, parvalues=None, parlimits=None,
            parlimited=None, parfixed=None, parerror=None, partied=None, fitunits=None):

        self.modelfunc = modelfunc
        self.npars = npars
        self.parnames = parnames
        self.fitunits = fitunits

        temp_pardict = dict([(varname, np.zeros(self.npars, dtype='bool')) if locals()[varname] is None else (varname, locals()[varname])
            for varname in str.split("parnames,parvalues,parlimits,parlimited,parfixed,parerror,partied",",")])

        # generate the parinfo dict
        # note that 'tied' must be a blank string (i.e. ""), not False, if it is not set
        # parlimited, parfixed, and parlimits are all two-element items (tuples or lists)
        self.parinfo = [ {'n':ii,
            'value':temp_pardict['parvalues'][ii],
            'limits':temp_pardict['parlimits'][ii],
            'limited':temp_pardict['parlimited'][ii],
            'fixed':temp_pardict['parfixed'][ii],
            'parname':temp_pardict['parnames'][ii],
            'error':temp_pardict['parerror'][ii],
            'tied':temp_pardict['partied'][ii] if temp_pardict['partied'][ii] else ""} 
            for ii in xrange(self.npars)]

    def mpfitfun(self,x,y,err):
        if err is None:
            def f(p,fjac=None): return [0,(y-self.modelfunc(x,*p))]
        else:
            def f(p,fjac=None): return [0,(y-self.modelfunc(x,*p))/err]
        return f

    def __call__(self, xax, data, err=None, guesses=[], quiet=True, shh=True, veryverbose=False, **kwargs):
        """
        Run the fitter
        """

        if len(guesses) == self.npars:
            for par,guess in zip(self.parinfo,guesses):
                par['value'] = guess

        if hasattr(xax,'convert_to_unit') and self.fitunits is not None:
            # some models will depend on the input units.  For these, pass in an X-axis in those units
            # (gaussian, voigt, lorentz profiles should not depend on units.  Ammonia, formaldehyde,
            # H-alpha, etc. should)
            xax = copy.copy(xax)
            xax.convert_to_unit(self.fitunits, quiet=quiet)

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


