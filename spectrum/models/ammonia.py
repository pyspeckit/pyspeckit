"""
Ammonia inversion transition TKIN fitter translated from Erik Rosolowsky's
http://svn.ok.ubc.ca/svn/signals/nh3fit/
"""
import numpy as np
from mpfit import mpfit
import spectrum.units as units
import matplotlib.cbook as mpcb

line_names = ['oneone','twotwo','threethree','fourfour']

freq_dict = { 
    'oneone':     23.694506e9,
    'twotwo':     23.722633335e9,
    'threethree': 23.8701296e9,
    'fourfour':   24.1394169e9,
    }
aval_dict = {
    'oneone':     1.712e-7,  #64*!pi**4/(3*h*c**3)*nu11**3*mu0**2*(1/2.)
    'twotwo':     2.291e-7,  #64*!pi**4/(3*h*c**3)*nu22**3*mu0**2*(2/3.)
    'threethree': 2.625e-7,  #64*!pi**4/(3*h*c**3)*nu33**3*mu0**2*(3/4.)
    'fourfour':   3.167e-7,  #64*!pi**4/(3*h*c**3)*nu44**3*mu0**2*(4/5.)
    }
ortho_dict = {
    'oneone':     False,
    'twotwo':     False,
    'threethree': True,
    'fourfour':   False,
    }

voff_lines_dict = {
    'oneone': [19.8513, 19.3159, 7.88669, 7.46967, 7.35132, 0.460409, 0.322042,
        -0.0751680, -0.213003, 0.311034, 0.192266, -0.132382, -0.250923, -7.23349,
        -7.37280, -7.81526, -19.4117, -19.5500],
    'twotwo':[26.5263, 26.0111, 25.9505, 16.3917, 16.3793, 15.8642, 0.562503,
        0.528408, 0.523745, 0.0132820, -0.00379100, -0.0132820, -0.501831,
        -0.531340, -0.589080, -15.8547, -16.3698, -16.3822, -25.9505, -26.0111,
        -26.5263],
    'threethree':[29.195098, 29.044147, 28.941877, 28.911408, 21.234827,
        21.214619, 21.136387, 21.087456, 1.005122, 0.806082, 0.778062,
        0.628569, 0.016754, -0.005589, -0.013401, -0.639734, -0.744554,
        -1.031924, -21.125222, -21.203441, -21.223649, -21.076291, -28.908067,
        -28.938523, -29.040794, -29.191744],
    'fourfour':[  0.        , -30.49783692,  30.49783692,   0., 24.25907811,
        -24.25907811,   0.        ]
                      }

tau_wts_dict = {
    'oneone': [0.0740740, 0.148148, 0.0925930, 0.166667, 0.0185190, 0.0370370,
        0.0185190, 0.0185190, 0.0925930, 0.0333330, 0.300000, 0.466667,
        0.0333330, 0.0925930, 0.0185190, 0.166667, 0.0740740, 0.148148],
    'twotwo': [0.00418600, 0.0376740, 0.0209300, 0.0372090, 0.0260470,
        0.00186000, 0.0209300, 0.0116280, 0.0106310, 0.267442, 0.499668,
        0.146512, 0.0116280, 0.0106310, 0.0209300, 0.00186000, 0.0260470,
        0.0372090, 0.0209300, 0.0376740, 0.00418600],
    'threethree': [0.012263, 0.008409, 0.003434, 0.005494, 0.006652, 0.008852,
        0.004967, 0.011589, 0.019228, 0.010387, 0.010820, 0.009482, 0.293302,
        0.459109, 0.177372, 0.009482, 0.010820, 0.019228, 0.004967, 0.008852,
        0.006652, 0.011589, 0.005494, 0.003434, 0.008409, 0.012263],
    'fourfour': [0.2431, 0.0162, 0.0162, 0.3008, 0.0163, 0.0163, 0.3911]}

class ammonia_model(object):

    def __init__(self):
        self.npeaks = 1
        self.npars = 6
        pass

    def ammonia(self, xarr, xunits='GHz', tkin=20, tex=20, Ntot=1e10, width=1,
            xoff_v=0.0, fortho=0.5, tau11=None, line='oneone'):
        """
        Generate a model Ammonia spectrum based on input temperatures, column, and
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

        if tex > tkin: # cannot have Tex > Tkin
            tex = tkin 

        ckms = 2.99792458e5
        ccms = ckms*1e5
        g1 = 1                
        g2 = 1                
        h = 6.6260693e-27     
        kb = 1.3806505e-16     
        mu0 = 1.476e-18               # Dipole Moment in cgs (1.476 Debeye)
      
        # Generate Partition Functions  
        nlevs = 51
        jv=np.arange(nlevs)
        ortho = jv % 3 == 0
        para = True-ortho
        Jpara = jv[para]
        Jortho = jv[ortho]
        Brot = 298117.06e6
        Crot = 186726.36e6
        Zpara = (2*Jpara+1)*np.exp(-h*(Brot*Jpara*(Jpara+1)+
            (Crot-Brot)*Jpara**2)/(kb*tkin))
        Zortho = 2*(2*Jortho+1)*np.exp(-h*(Brot*Jortho*(Jortho+1)+
            (Crot-Brot)*Jortho**2)/(kb*tkin))

        runspec = np.zeros(len(xarr))
        
        tau_dict = {}
        for ii,linename in enumerate(line_names):
            orthoparafrac = fortho if ortho_dict[linename] else (1-fortho)
            Z = Zortho if ortho_dict[linename] else Zpara
            tau_dict[linename] = (Ntot * orthoparafrac * Z[ii]/(Z.sum()) / ( 1
                + np.exp(-h*freq_dict[linename]/(kb*tkin) )) * ccms**2 /
                (8*np.pi*freq_dict[linename]**2) * aval_dict[linename]*
                (1-np.exp(-h*freq_dict[linename]/(kb*tex))) /
                (width/ckms*freq_dict[linename]*np.sqrt(2*np.pi)) )
      
            voff_lines = np.array(voff_lines_dict[linename])
            tau_wts = np.array(tau_wts_dict[linename])
      
            lines = (1-voff_lines/ckms)*freq_dict[linename]/1e9
            tau_wts = tau_wts / (tau_wts).sum()
            nuwidth = np.abs(width/ckms*lines)
            nuoff = xoff_v/ckms*lines
      
            # tau array
            tauprof = np.zeros(len(xarr))
            for kk,no in enumerate(nuoff):
              tauprof += tau_dict[linename]*tau_wts[kk]*\
                        np.exp(-(xarr+no-lines[kk])**2/(2*nuwidth[kk]**2))
      
            T0 = (h*xarr*1e9/kb)
            runspec = (T0/(np.exp(T0/tex)-1)-T0/(np.exp(T0/2.73)-1))*(1-np.exp(-tauprof))+runspec
            if runspec.min() < 0:
                raise ValueError("Model dropped below zero.  That is not possible normally.")
      
        return runspec

    def n_ammonia(self, pars=None,**kwargs):
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

        def L(x):
            v = np.zeros(len(x))
            for i in range(len(tkin)):
                v += self.ammonia(x,tkin=tkin[i],tex=tex[i],Ntot=Ntot[i],width=width[i],xoff_v=xoff_v[i],fortho=fortho[i],**kwargs)
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
            if err == None:
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
