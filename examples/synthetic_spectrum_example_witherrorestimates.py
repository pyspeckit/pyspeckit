import numpy as np
from astropy import units as u
import itertools
import pyspeckit
import scipy.stats

from astropy.utils.console import ProgressBar

import pylab as pl
pl.close('all')
pl.figure(1).clf()

xaxis = np.linspace(-50.,150.,100)
sigma = 10.
center = 50.
synth_data = np.exp(-(xaxis-center)**2/(sigma**2 * 2.))

# Add noise
stddev = 0.1
noise = np.random.randn(xaxis.size)*stddev
error = stddev*np.ones_like(synth_data)
data = noise+synth_data

# this will give a "blank header" warning, which is fine
sp = pyspeckit.Spectrum(data=data, error=error, xarr=xaxis,
                        xarrkwargs={'unit':'km/s'},
                        unit=u.erg/u.s/u.cm**2/u.AA,
                        header={},
                       )

sp.plotter(figure=pl.figure(1))

# Fit with automatic guesses
sp.specfit(fittype='gaussian')

# Fit with input guesses
# The guesses initialize the fitter
# This approach uses the 0th, 1st, and 2nd moments
amplitude_guess = data.max()
center_guess = (data*xaxis).sum()/data.sum()
width_guess = (data.sum() / amplitude_guess / (2*np.pi))**0.5
guesses = [amplitude_guess, center_guess, width_guess]
sp.specfit(fittype='gaussian', guesses=guesses)

sp.plotter(errstyle='fill')
sp.specfit.plot_fit()


sp.plotter.savefig('oned_gaussfit_example.pdf')


# do a deeper examination of the likelihood function

def chi2(sp, pars):
    """
    Given a spectrum and some parameters, calculate the chi^2 value
    """
    return ((sp.specfit.get_model_frompars(sp.xarr, pars) -
             sp.specfit.spectofit)**2 /
            (sp.specfit.errspec**2)
           ).sum()

def replace(lst, eltid, value):
    """
    Helper function to make new parameter value lists below
    """
    lstcopy = list(lst)
    lstcopy[eltid] = value
    return lstcopy

nsteps = {'AMPLITUDE0':21, 'SHIFT0':19, 'WIDTH0': 20}
par_edges = [np.linspace(par.value-par.error*2,
                         par.value+par.error*2,
                         nsteps[par.parname])
             for par in sp.specfit.parinfo]
par_cube = np.meshgrid(*par_edges, indexing='ij')
params = np.array(list(zip(map(np.ravel, par_cube)))).squeeze().T

par_like_cube = np.array([chi2(sp, pars) for pars in
                          ProgressBar(params)]).reshape(par_cube[0].shape)

fig = pl.figure(2, figsize=(12,12))
fig.clf()

for ii,par in enumerate(sp.specfit.parinfo):
    other_axes = tuple([x for x in (0,1,2) if x != ii])
    par_likes = par_like_cube.min(axis=other_axes)
    par_vals = par_edges[ii]

    pct68 = scipy.stats.chi2.cdf(1, 1)

    # this is for the full 3-parameter case
    delta_chi2 = scipy.stats.chi2.ppf(pct68, sp.specfit.fitter.npars)

    # this is the one we're actually using (which is just 1) because we're
    # marginalizing over the other parameters
    delta_chi2 = scipy.stats.chi2.ppf(pct68, 1)

    ax = fig.add_subplot(3,3,1+ii*4)
    ax.plot(par_vals, par_likes)
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.hlines(par_likes.min() + delta_chi2, xmin, xmax, 'k', '--')
    ax.vlines([par.value - par.error,
               par.value + par.error,], ymin, ymax, 'k', '--')
    ax.set_ylabel("$\Delta \chi^2$")
    ax.set_xlabel("{0} Value".format(par.parname))
    pl.setp(ax.get_xticklabels(), rotation=90)

    ax.set_xlim(xmin, xmax,)
    ax.set_ylim(ymin, ymax,)

fig.tight_layout()


plotinds = {0:4, 1:7, 2:8}
for ii,((ind1,ind2),(par1,par2)) in enumerate(zip(itertools.combinations([0,1,2],2),
                                                  itertools.combinations(sp.specfit.parinfo,2)
                                                 )):

    # list(itertools.combinations((1,2,3),2)) = [(1, 2), (1, 3), (2, 3)]
    other_axis, = {0,1,2} - {ind1,ind2}

    par_likes = par_like_cube.min(axis=other_axis)
    p1vals = par_edges[ind1]
    p2vals = par_edges[ind2]

    pct68 = 1-scipy.stats.norm.sf(1)*2
    pct95 = 1-scipy.stats.norm.sf(2)*2
    pct997 = 1-scipy.stats.norm.sf(3)*2

    # marginalize over only one parameter; we now have two free parameters
    delta_chi2_68 = scipy.stats.chi2.ppf(pct68, 2)
    delta_chi2_95 = scipy.stats.chi2.ppf(pct95, 2)
    delta_chi2_997 = scipy.stats.chi2.ppf(pct997, 2)

    ax = fig.add_subplot(3,3,plotinds[ii])
    ax.contour(p1vals, p2vals, par_likes.T,
               levels=[par_likes.min()+1,
                       par_likes.min()+delta_chi2_68,
                       par_likes.min()+delta_chi2_95,
                       par_likes.min()+delta_chi2_997],
               colors=['k', 'r', 'b', 'g'],
               linestyles=['--', '-', ':', '-.'],
              )
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.plot(par1.value, par2.value, 'kx')
    ax.hlines([par2.value - par2.error,
               par2.value + par2.error,], xmin, xmax, 'k', '--')
    ax.vlines([par1.value - par1.error,
               par1.value + par1.error,], ymin, ymax, 'k', '--')
    ax.set_xlabel("{0} Value".format(par1.parname))
    ax.set_ylabel("{0} Value".format(par2.parname))
    pl.setp(ax.get_xticklabels(), rotation=90)

    ax.set_xlim(xmin, xmax,)
    ax.set_ylim(ymin, ymax,)

fig.tight_layout()

fig.savefig("error_estimate_demonstration.pdf")
