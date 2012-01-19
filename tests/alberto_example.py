#! /usr/bin/env python
import pyspeckit
from pylab import *

if not 'savedir' in globals():
    savedir = ''

c = 2.99792458e18 # A/s

direc = './'
files = ['alberto_example.txt'] # ready for multiple objects

X = []
spec = []

#def f(x):
#   return x

for i in range(len(files)):
    # reading every object in files and converting to my favorite units
    z0 = 0.985
    X.append(loadtxt(direc+files[i]))
    wav= X[i][:,0] # Angstrom
    flux = X[i][:,0]**2/c*1e29*X[i][:,1] # microJy
    fluxerr = X[i][:,0]**2/c*1e29*X[i][:,2] # microJy

    spec.append(pyspeckit.Spectrum(xarr=wav,data=flux,err=fluxerr,units='Jy',xarrkwargs={'unit':'angstroms'}))
    spec[i].units = 'Jy'
    spec[i].xarr.frame = 'obs'
    spec[i].xarr.redshift = z0
    #my_fitter = f.SpectralModel(hill5_model,3,parnames=['amplitude','wav','width'],parlimited=[(False,False),(False,False),(False,False)],parlimits=[(0,0),(0,0),(0,0),(0,0),(0,0)],shortvarnames=('A','\lambda','\sigma'))

    # prepare spectrum and remove overall continuum
    spec[i].plotter()
    spec[i].baseline(subtract=False)
    spec[i].baseline()
    #spec[i].baseline.selectregion(xmin=9000, xmax=10500, xtype='wcs', highlight=True)
    spec[i].plotter.refresh()
    spec[i].plotter.figure.savefig(savedir+'example1a.png')
    # remove continuum around the emission lines
    spec[i].crop(9000,10500)
    spec[i].baseline(subtract=False)
    spec[i].baseline()
    spec[i].plotter.refresh()
    spec[i].plotter.figure.savefig(savedir+'example1b.png')

    # define my emission lines
    Hbeta = spec[i].speclines.optical.lines['H_beta'][0]
    OIIIa = spec[i].speclines.optical.lines['OIIIa'][0]
    OIIIb = spec[i].speclines.optical.lines['OIIIb'][0]

    # [Amplitude, wavelength, width]
    guesses = [10,Hbeta*(1+z0),5,30,OIIIa*(1+z0),50,90,OIIIb*(1+z0),50]
    # fitting
    spec[i].specfit(guesses=guesses,quiet=False,annotate=False)

    ## one way 'tied' can be implemented:
    #parinfo = spec[i].specfit.parinfo
    #parinfo[7]['tied'] = 'p[4]+%f' % (guesses[7]-guesses[4])
    #parinfo[1]['tied'] = 'p[4]-%f' % (guesses[4]-guesses[1])
    #spec[i].specfit(guesses=guesses,quiet=False,annotate=False,parinfo=parinfo)

    # here's another way:
    parinfo = spec[i].specfit.parinfo
    tied = [''] * 10 # make an empty 'tied' list
    tied[7] = 'p[1]+%f' % (guesses[7]-guesses[4])
    tied[4] = 'p[1]+%f' % (guesses[4]-guesses[1])
    spec[i].specfit(guesses=guesses,quiet=False,annotate=False,tied=tied)

    # you probably want to do this - iteratively determine continuum/baseline
    spec[i].baseline(excludefit=True) # fit the baseline excluding the emission lines
    # then refit
    spec[i].specfit(guesses=guesses,quiet=False,annotate=False,parinfo=parinfo)
    # cross-checking redshift
    OIIIb_obs = spec[i].specfit.modelpars[-2]
    print 'Our guess for the redshift was z = %g' % z0
    print 'The redshift, as derived by the line shift, is z = %g' % ((OIIIb_obs/OIIIb)-1)

    spec[i].specfit.shift_pars('rest')
    
    # measuring line properties
    spec[i].measure(z=z0)

    # Make some nice axis labels
    xticks(fontsize=18)
    yticks(fontsize=18)
    spec[i].plotter.axis.set_xlabel(r'$\lambda_{obs}$ [$\AA$]',fontsize=18)
    spec[i].plotter.axis.set_ylabel(r'$F_{\nu}$ [$\mu$Jy]',fontsize=18)
    spec[i].plotter.refresh()
    spec[i].plotter.figure.savefig(savedir+'example1c.png')

