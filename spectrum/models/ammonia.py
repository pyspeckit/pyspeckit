import numpy as np

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

def ammonia(xarr, xtype='frequency', xunits='GHz', tkin=20, tex=20, Ntot=1e10, width=1,
        xcen=0.0, fortho=0.5):
    """
    Generate a model Ammonia spectrum based on input temperatures, column, and
    gaussian parameters
    """

    if xunits == 'Hz':
        xarr /= 1e9
    if xtype != 'frequency':
        raise Exception( "Error: convert to frequency first" )

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
    for linename in line_names:
        orthoparafrac = fortho if ortho_dict[linename] else (1-fortho)
        tau_dict[linename] = (Ntot * orthoparafrac * Zpara[0]/(Zpara.sum()) / ( 1
            + np.exp(-h*freq_dict[linename]/(kb*tkin) )) * ccms**2 /
            (8*np.pi*freq_dict[linename]**2) * aval_dict[linename]*
            (1-np.exp(-h*freq_dict[linename]/(kb*tex))) /
            (width/ckms*freq_dict[linename]*np.sqrt(2*np.pi)) )
  
        voff_lines = np.array(voff_lines_dict['oneone'])
        tau_wts = np.array(tau_wts_dict['oneone'])
  
        lines = (1-voff_lines/ckms)*freq_dict[linename]/1e9
        tau_wts = tau_wts / (tau_wts).sum()
        nuwidth = np.abs(width/ckms*lines)
        nuoff = xcen/ckms*lines
  
        # tau array
        tauprof = np.zeros(len(xarr))
        for kk,no in enumerate(nuoff):
          tauprof += tau_dict[linename]*tau_wts[kk]*\
                    np.exp(-(xarr-no-lines[kk])**2/(2*nuwidth[kk]**2))
  
        T0 = (h*xarr*1e9/kb)
        runspec = (T0/(np.exp(T0/tex)-1)-T0/(np.exp(T0/2.73)-1))*(1-np.exp(-tauprof))+runspec
  
    return runspec

def n_ammonia(pars=None,**kwargs):
    """
    Returns a function that sums over N ammonia line profiles, where N is the length of
    tkin,tex,Ntot,width,xcen,fortho *OR* N = len(pars) / 6

    The background "height" is assumed to be zero (you must "baseline" your
    spectrum before fitting)

    pars  - a list with len(pars) = 6n, assuming tkin,tex,Ntot,width,xcen,fortho repeated
    """
    if len(pars) % 6 == 0:
        tkin = [pars[ii] for ii in xrange(0,len(pars),6)]
        tex = [pars[ii] for ii in xrange(1,len(pars),6)]
        Ntot = [pars[ii] for ii in xrange(2,len(pars),6)]
        width = [pars[ii] for ii in xrange(3,len(pars),6)]
        xcen = [pars[ii] for ii in xrange(4,len(pars),6)]
        fortho = [pars[ii] for ii in xrange(5,len(pars),6)]
    elif not(len(dx) == len(width) == len(a)):
        raise ValueError("Wrong array lengths! dx: %i  width %i  a: %i" % (len(dx),len(width),len(a)))

    def L(x):
        v = numpy.zeros(len(x))
        for i in range(len(dx)):
            v += ammonia(x,tkin=tkin[i],tex=tex[i],Ntot=Ntot[i],width=width[i],xcen=xcen[i],fortho=fortho[i],**kwargs)
        return v
    return L

def multinh3fit(xax, data, nnh3=1, err=None, params=[20,20,1e10,1.0,0.0,0.5],
        fixed=[False,False,False,False,False,False],
        limitedmin=[True,True,True,True,False,True],
        limitedmax=[False,False,False,False,False,True], minpars=[2.73,0,0,0,0,0],
        maxpars=[0,0,0,0,0,1], quiet=True, shh=True, veryverbose=False):
    """
    Fit multiple nh3 profiles

    Inputs:
       xax - x axis
       data - y axis
       nnh3 - How many nh3 profiles to fit?  Default 1 (this could supersede onedgaussfit)
       err - error corresponding to data

     These parameters need to have length = 6*nnh3.  If nnh3 > 1 and length = 6, they will
     be replicated nnh3 times, otherwise they will be reset to defaults:
       params - Fit parameters: [amplitude, offset, Gfwhm, Lfwhm] * nnh3
              If len(params) % 6 == 0, nnh3 will be set to len(params) / 6
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

    npars = 6

    if len(params) != nnh3 and (len(params) / npars) > nnh3:
        nnh3 = len(params) / npars 

    if isinstance(params,numpy.ndarray): params=params.tolist()

    # make sure all various things are the right length; if they're not, fix them using the defaults
    for parlist in (params,fixed,limitedmin,limitedmax,minpars,maxpars):
        if len(parlist) != npars*nnh3:
            # if you leave the defaults, or enter something that can be multiplied by 3 to get to the
            # right number of gaussians, it will just replicate
            if len(parlist) == npars: 
                parlist *= nnh3 
            elif parlist==params:
                parlist[:] = [20,20,1e10,1.0,0.0,0.5] * nnh3
            elif parlist==fixed:
                parlist[:] = [False,False,False,False,False,False] * nnh3
            elif parlist==limitedmax:
                parlist[:] = [False,False,False,False,False,True] * nnh3
            elif parlist==limitedmin:
                parlist[:] = [True,True,True,True,False,True] * nnh3
            elif parlist==minpars:
                parlist[:] = [2.73,0,0,0,0,0] * nnh3
            elif parlist==maxpars:
                parlist[:] = [0,0,0,0,0,1] * nnh3

    def mpfitfun(x,y,err):
        if err == None:
            def f(p,fjac=None): return [0,(y-n_ammonia(pars=p)(x))]
        else:
            def f(p,fjac=None): return [0,(y-n_ammonia(pars=p)(x))/err]
        return f

    parnames = {0:"TKIN",1:"TEX",2:"NTOT",3:"WIDTH",4:"XCEN",5:"FORTHO"}

    parinfo = [ {'n':ii, 'value':params[ii],
        'limits':[minpars[ii],maxpars[ii]],
        'limited':[limitedmin[ii],limitedmax[ii]], 'fixed':fixed[ii],
        'parname':parnames[ii%npars]+str(ii%npars), 'error':ii} 
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

    return mpp,n_ammonia(pars=mpp)(xax),mpperr,chi2

