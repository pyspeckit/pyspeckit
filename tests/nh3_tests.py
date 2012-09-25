import pyspeckit
from pyspeckit.spectrum.models import ammonia
import time
import numpy as np
import itertools
import pylab
try:
    import astropy.io.fits as pyfits
    import astropy.wcs as pywcs
except ImportError:
    import pyfits
    import pywcs
try: 
    import progressbar
    widgets = [progressbar.FormatLabel('Processed: %(value)d offsets in %(elapsed)s)'), progressbar.Percentage()]
    progress = progressbar.ProgressBar(widgets=widgets)
except ImportError:
    def progress(x):
        yield x

xarr11 = pyspeckit.units.SpectroscopicAxis(np.arange(-30,30,0.4),
        units='km/s', refX=23.694496, refX_units='GHz', frame='LSRK',
        xtype='Frequency')

xarr22 = pyspeckit.units.SpectroscopicAxis(np.arange(-30,30,0.4),
        units='km/s', refX=23.722633335, refX_units='GHz', frame='LSRK',
        xtype='Frequency')

xarr = pyspeckit.units.SpectroscopicAxes([xarr11.as_unit('GHz'),xarr22.as_unit('GHz')])

test_parameters = {
        'TKIN': [7.5,10,12.5,15,17.5,20,25,30,35,40],
        'Tex_frac': np.linspace(0.1,1,4),
        'ntot': np.linspace(14,17,5),
        'width': [0.5],
        }

parnames = ['TKIN','Tex_frac','ntot','width']

parlist = (list(itertools.product(*[list(test_parameters[k]) for k in parnames])))

# typical error is ~0.4 K in 0.4 km/s channels [ERRORS ARE OVERESTIMATED IN FITTER...]

noise = 0.4
err = np.ones(300,dtype='float') * noise 

# CASA-like convenience definition
F = False
T = True

ins = []
outs = []
errors = []
taus = []
noises = []
signoise = []

if True:
    tk = 15
    tx = 7.5
    nt = 15
    sig = 0.5
np.seterr(all='ignore')
for tk in test_parameters['TKIN']:
    for txf in test_parameters['Tex_frac']:
        for nt in test_parameters['ntot']:
            for sig in test_parameters['width']:
#for tk, txf, nt, sig in (parlist):

                t0 = time.time()
                tx = max([txf*tk,3]) # tx can't be too close to T(CMB)

                ins.append([tk,tx,nt,sig])

                fake_spec = ammonia.ammonia(xarr, tkin=tk, tex=tx, ntot=nt, width=sig, xoff_v=0,
                        fortho=0.5)
                real_tau = ammonia.ammonia(xarr, tkin=tk, tex=tx, ntot=nt, width=sig, xoff_v=0,
                        fortho=0.5, return_tau=True)

                if fake_spec.max()*3 < noise:
                    noises.append(fake_spec.max()/3.)
                    fake_spec += np.random.randn(fake_spec.size) * fake_spec.max()/3.
                    signoise.append(3)
                else:
                    noises.append(noise)
                    signoise.append(fake_spec.max()/noise)
                    fake_spec += np.random.randn(fake_spec.size) * noise

                sp = pyspeckit.Spectrum(xarr=xarr, data=fake_spec, err=err,units='K', header=pyfits.Header({'BUNIT':'K'}))
                
                # sp.plotter()

                sp.specfit(fittype='ammonia', multifit=True, guesses=[tk,tx,nt,sig,0,0.5],
                        fixed=[F,F,F,F,F,T], 
                        limitedmax=[T,T,T,T,T,T],
                        maxpars=[100,100,25,5,20,1],
                        limitedmin=[T,T,T,T,T,T],
                        minpars=[2.73,2.73,10,0,-20,0],
                        use_lmfit=True,
                        verbose=False)

                outs.append(sp.specfit.parinfo.values[:4])
                errors.append(sp.specfit.parinfo.errors[:4])
                print " ".join(["%7s" % round(x,1) for x in ins[-1]]),"|",
                print " ".join(["%7s" % round(x,1) for x in outs[-1]]),"|",
                print " ".join(["%7s" % round((x-y),1) for x,y in zip(ins[-1],outs[-1])]),"|",
                print " ".join(["%7s" % round((x),1) for x in errors[-1]]),"|",

                fit_tau = ammonia.ammonia(xarr, tkin=outs[-1][0], tex=outs[-1][1], ntot=outs[-1][2], width=outs[-1][3], xoff_v=0,
                        fortho=0.5, return_tau=True)

                taus.append([real_tau['oneone'],real_tau['twotwo'],
                        fit_tau['oneone'],fit_tau['twotwo']])
                print " ".join(['%7s' % round(x,1) for x in taus[-1]]),"|",
                print "%7.1f" % (time.time()-t0)

    insarr = np.array(ins)
    outsarr = np.array(outs)
    tausarr = np.array(taus)
    errarr = np.array(errors)
    signarr = np.array(signoise)
    pylab.figure()
    for te in np.unique(insarr[:,1]):
        okT = insarr[:,0] == tk
        texp5 = insarr[:,1] == te
        widthp5 = insarr[:,3] == 0.5
        OK = okT*texp5*widthp5
        if OK.sum() == 0:
            continue
        pylab.plot(insarr[OK,2], insarr[OK,0], 'o', label='Real T$_K$ (T$_X$=%0.1f)' % te)
        pylab.errorbar(outsarr[OK,2], outsarr[OK,0], yerr=errarr[OK,0], xerr=errarr[OK,2], linestyle='none', label='Fitted T$_K$ (T$_X$=%0.1f)' % te)
    pylab.axis([13.5,17.5,tk*0.5,tk*1.5])
    pylab.xlabel("Log Column Density")
    pylab.ylabel("Kinetic Temperature")
    pylab.title("T=%0.1f" % tk)
    pylab.legend(loc='best')

    pylab.figure()
    for te in np.unique(insarr[:,1]):
        okT = insarr[:,0] == tk
        texp5 = insarr[:,1] == te
        widthp5 = insarr[:,3] == 0.5
        OK = okT*texp5*widthp5
        if OK.sum() == 0:
            continue
        pylab.plot(signarr[OK], insarr[OK,0], 'o', label='Real T$_K$ (T$_X$=%0.1f)' % te)
        pylab.errorbar(signarr[OK], outsarr[OK,0], yerr=errarr[OK,0], linestyle='none', label='Fitted T$_K$ (T$_X$=%0.1f)' % te)
    pylab.axis([3,100,tk*0.5,tk*1.5])
    pylab.xlabel("Signal to Noise Ratio")
    pylab.ylabel("Kinetic Temperature")
    pylab.title("T=%0.1f" % tk)
    pylab.legend(loc='best')

    pylab.figure()
    for te in np.unique(insarr[:,1]):
        okT = insarr[:,0] == tk
        texp5 = insarr[:,1] == te
        widthp5 = insarr[:,3] == 0.5
        OK = okT*texp5*widthp5
        if OK.sum() == 0:
            continue
        pylab.scatter(signarr[OK], insarr[OK,2], s=tausarr[OK,0], marker='o', label='Real T$_K$ (T$_X$=%0.1f)' % te)
        pylab.errorbar(signarr[OK], outsarr[OK,2], yerr=errarr[OK,2], linestyle='none', label='Fitted T$_K$ (T$_X$=%0.1f)' % te)
    pylab.axis([3,100,13.5,17.5])
    pylab.xlabel("Signal to Noise Ratio")
    pylab.ylabel("Column Density")
    pylab.title("T=%0.1f" % tk)
    pylab.legend(loc='best')

