import numpy as np
from astropy import units as u
import itertools
from operator import itemgetter
import pyspeckit
import scipy.stats

from astropy.utils.console import ProgressBar

from pyspeckit.spectrum.models import ammonia, ammonia_constants

import pylab as pl
pl.close('all')
pl.figure(1).clf()

oneonefreq = ammonia_constants.freq_dict['oneone']
twotwofreq = ammonia_constants.freq_dict['twotwo']
# create an axis that covers the 1-1 and 2-2 inversion lines
xaxis1 = pyspeckit.units.SpectroscopicAxis(np.linspace(oneonefreq*(1-50/3e5),
                                                       oneonefreq*(1+50/3e5),
                                                       80.), unit=u.Hz,
                                           velocity_convention='radio',
                                           refX=oneonefreq*u.Hz)
xaxis2 = pyspeckit.units.SpectroscopicAxis(np.linspace(twotwofreq*(1-50/3e5),
                                                       twotwofreq*(1+50/3e5),
                                                       80.), unit=u.Hz,
                                           velocity_convention='radio',
                                           refX=twotwofreq*u.Hz)

sigma = 2.
center = 0.
trot = 35.
ntot = 15.5
tex = 35. # Adopting an LTE model
synth_data1 = ammonia.ammonia(xaxis1, trot=trot, tex=tex, width=sigma,
                              xoff_v=center, ntot=ntot)
synth_data2 = ammonia.ammonia(xaxis2, trot=trot, tex=tex, width=sigma,
                              xoff_v=center, ntot=ntot)

# Add noise
stddev = 0.4
noise = np.random.randn(xaxis1.size)*stddev
error = stddev*np.ones_like(synth_data1)
data1 = noise+synth_data1
noise = np.random.randn(xaxis2.size)*stddev
data2 = noise+synth_data2

# this will give a "blank header" warning, which is fine
sp1 = pyspeckit.Spectrum(data=data1, error=error, xarr=xaxis1,
                         xarrkwargs={'unit':'km/s', 'refX': oneonefreq*u.Hz,},
                         unit=u.K)
sp2 = pyspeckit.Spectrum(data=data2, error=error, xarr=xaxis2,
                         xarrkwargs={'unit':'km/s', 'refX': twotwofreq*u.Hz,},
                         unit=u.K)
sp = pyspeckit.Spectra([sp1,sp2])

#sp.plotter(figure=pl.figure(1), errstyle='fill')

# fit with some vague initial guesses and a fixed ortho/para
# fraction of 1 (it is unconstrained in our data)
sp.specfit(fittype='ammonia',
           #guesses=[20, 15, 15.5, 3, 2, 1],
           guesses=[30, 25, 15.5, 3, 2, 0],
           fixed=[False,False,False,False,False,True])

# do it again to set the errors neatly
sp.specfit(fittype='ammonia',
           #guesses=[20, 15, 15.5, 3, 2, 1],
           guesses=sp.specfit.parinfo.values,
           fixed=[False,False,False,False,False,True])


pyspeckit.wrappers.fitnh3.plot_nh3({'oneone':sp1, 'twotwo': sp2}, sp,
                                  show_hyperfine_components=False,
                                  errstyle='fill')
#sp.plotter(errstyle='fill')
#sp.specfit.plot_fit()

sp.plotter.figure.subplots_adjust(hspace=0.5)

sp.plotter.savefig('oned_ammonia_LTE_fit_example.pdf')

sp.xarr.convert_to_unit(u.GHz)

# do a deeper examination of the likelihood function

def chi2(sp, pars):
    """
    Given a spectrum and some parameters, calculate the chi^2 value
    """
    pars = list(pars) + [0]
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

nsteps = {'width0': 7, 'tex0': 21, 'trot0': 20, 'ntot0': 19}
par_edges = [np.linspace(par.value-par.error*2,
                         par.value+par.error*2,
                         nsteps[par.parname])
             for par in sp.specfit.parinfo
             if not par.fixed and not par.tied
             and not par.parname == 'xoff_v0'
            ]
par_cube = np.meshgrid(*par_edges, indexing='ij')
params = np.array(list(zip(map(np.ravel, par_cube)))).squeeze().T

xoff_v0 = sp.specfit.parinfo['xoff_v0'].value

if 'par_like_cube' not in locals():
    # allow re-run of plot code without re-rerunning the chi^2 calculations
    par_like_cube = np.array([chi2(sp, pars.tolist()+[xoff_v0]) for pars in
                              ProgressBar(params)]).reshape(par_cube[0].shape)


fig = pl.figure(2, figsize=(8,8))
fig.clf()
fontsize = 8

for ii,par in enumerate(sp.specfit.parinfo):
    if par.fixed or par.parname == 'xoff_v0':
        continue

    other_axes = tuple([x for x in (0,1,2,3) if x != ii])
    par_likes = par_like_cube.min(axis=other_axes)
    par_vals = par_edges[ii]

    pct68 = scipy.stats.chi2.cdf(1, 1)

    # this is for the full 3-parameter case
    delta_chi2 = scipy.stats.chi2.ppf(pct68, sp.specfit.fitter.npars - sum(sp.specfit.parinfo.fixed))

    # this is the one we're actually using (which is just 1) because we're
    # marginalizing over the other parameters
    delta_chi2 = scipy.stats.chi2.ppf(pct68, 1)

    ax = fig.add_subplot(4,4,1+ii*5)
    ax.plot(par_vals, par_likes)
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.hlines(par_likes.min() + delta_chi2, xmin, xmax, 'k', '--')
    ax.vlines([par.value - par.error,
               par.value + par.error,], ymin, ymax, 'k', '--')
    ax.set_ylabel("$\Delta \chi^2$", fontsize=fontsize)
    ax.set_xlabel("{0} Value".format(par.parname), fontsize=fontsize)
    pl.setp(ax.get_xticklabels(), rotation=90, fontsize=fontsize)
    pl.setp(ax.get_yticklabels(), fontsize=fontsize)

    ax.set_xlim(xmin, xmax,)
    ax.set_ylim(ymin, ymax,)

fig.tight_layout()


nsigma_range = 2
nsteps = 15.

plotinds = {0:6, 1:11, 2:12, 3:16, 4:17, 5:18, 6:21, 7:22, 8:23, 9:24}
plotinds = {(0,1): 6, (0,2): 11, (1,2): 12, (0,3): 16, (1,3): 17, (2,3): 18,
            (0,4): 21, (1,4): 22, (2,4): 23, (3,4): 24, }
plotinds = {(0,1): 5, (0,2): 9, (1,2): 10, (0,3): 13, (1,3): 14, (2,3): 15}

for ii,((ind1,ind2),(par1,par2)) in enumerate(zip(itertools.combinations([0,1,2,3], 2),
                                                  itertools.combinations(sp.specfit.parinfo[:4], 2)
                                                 )):

    other_axis = tuple({0,1,2,3} - {ind1,ind2})

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

    ax = fig.add_subplot(4,4,plotinds[(ind1,ind2)])
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
    ax.set_xlabel("{0} Value".format(par1.parname), fontsize=fontsize)
    ax.set_ylabel("{0} Value".format(par2.parname), fontsize=fontsize)
    pl.setp(ax.get_xticklabels(), rotation=90, fontsize=fontsize)
    pl.setp(ax.get_yticklabels(), fontsize=fontsize)

    ax.set_xlim(xmin, xmax,)
    ax.set_ylim(ymin, ymax,)

fig.tight_layout()
fig.subplots_adjust(hspace=0.7, wspace=0.7)

fig.savefig("ammonia_LTE_default_error_estimate_demonstration.pdf")


# fit with restricted delta-T instead of restricted T

sp.specfit.Registry.add_fitter('ammonia_dtex',
                               ammonia.ammonia_model_restricted_tex(), 7,
                               multisingle='multi')

sp.specfit(fittype='ammonia_dtex',
           guesses=[30, 25, 15.5, 3, 2, 0, 5],
           fixed=[False,False,False,False,False,True,False])

sp.plotter.savefig('oned_ammonia_LTE_fit_example_deltaT.pdf')


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


nsteps = {'width0': 7, 'delta0': 21, 'trot0': 20, 'xoff_v0': 1,  'ntot0': 19}
par_edges = [np.linspace(0, 10/nsteps[par.parname], nsteps[par.parname])
             if par.parname == 'delta0' else
             np.linspace(par.value-par.error*2,
                         par.value+par.error*2,
                         nsteps[par.parname])
             for par in sp.specfit.parinfo
             if not par.fixed and not par.tied
             and par.parname != 'xoff_v0'
            ]
par_cube = np.meshgrid(*par_edges, indexing='ij')
assert par_cube[0].shape == (nsteps['trot0'], nsteps['ntot0'],
                             nsteps['width0'],
                             #nsteps['xoff_v0'],
                             nsteps['delta0'])
params = list(map(list, zip(*map(list, map(np.ravel, par_cube)))))

xoff_v0 = sp.specfit.parinfo['xoff_v0'].value

if 'par_like_cube_dtex' not in locals():
    # allow re-run of plot code without re-rerunning the chi^2 calculations
    par_like_cube_dtex = np.array([chi2(sp,
                                   pars[0:1] + [pars[0]+pars[-1]] + pars[1:4] +
                                   [0, xoff_v0] + pars[5:])
                              for pars in
                              ProgressBar(params)]).reshape(par_cube[0].shape)


pl.close(2)
# this will fail in ugly ways if run on a small screen (everything will be squished)
fig = pl.figure(2, figsize=(8,8))
fontsize = 8
fig.clf()

dd = 0
for ii,par in enumerate(sp.specfit.parinfo):
    if par.fixed or par.tied or par.parname == 'xoff_v0':
        continue

    other_axes = tuple([x for x in (0,1,2,3) if x != dd])
    par_likes = par_like_cube_dtex.min(axis=other_axes)
    par_vals = par_edges[dd]

    pct68 = scipy.stats.chi2.cdf(1, 1)

    # this is for the full 3-parameter case
    delta_chi2 = scipy.stats.chi2.ppf(pct68, sp.specfit.fitter.npars - sum(sp.specfit.parinfo.fixed))

    # this is the one we're actually using (which is just 1) because we're
    # marginalizing over the other parameters
    delta_chi2 = scipy.stats.chi2.ppf(pct68, 1)

    ax = fig.add_subplot(4,4,1+(dd)*5)
    ax.plot(par_vals, par_likes)
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.hlines(par_likes.min() + delta_chi2, xmin, xmax, 'k', '--')
    ax.vlines([par.value - par.error,
               par.value + par.error,], ymin, ymax, 'k', '--')
    ax.set_ylabel("$\Delta \chi^2$", fontsize=fontsize)
    ax.set_xlabel("{0} Value".format(par.parname), fontsize=fontsize)
    pl.setp(ax.get_xticklabels(), rotation=90, fontsize=fontsize)
    pl.setp(ax.get_yticklabels(), fontsize=fontsize)

    ax.set_xlim(xmin, xmax,)
    ax.set_ylim(ymin, ymax,)

    dd += 1

fig.tight_layout()


nsigma_range = 2
nsteps = 15.


plotinds = {0:6, 1:11, 2:12, 3:16, 4:17, 5:18, 6:21, 7:22, 8:23, 9:24}
plotinds = {(0,1): 6, (0,2): 11, (1,2): 12, (0,3): 16, (1,3): 17, (2,3): 18,
            (0,4): 21, (1,4): 22, (2,4): 23, (3,4): 24, }
plotinds = {(0,1): 5, (0,2): 9, (1,2): 10, (0,3): 13, (1,3): 14, (2,3): 15}
for ii,((ind1,ind2),(par1,par2)) in enumerate(zip(itertools.combinations([0,1,2,3], 2),
                                                  itertools.combinations(itemgetter(0,2,3,6)(sp.specfit.parinfo),2)
                                                 )):

    other_axis = tuple({0,1,2,3} - {ind1,ind2})

    par_likes = par_like_cube_dtex.min(axis=other_axis)
    p1vals = par_edges[ind1]
    p2vals = par_edges[ind2]


    pct68 = 1-scipy.stats.norm.sf(1)*2
    pct95 = 1-scipy.stats.norm.sf(2)*2
    pct997 = 1-scipy.stats.norm.sf(3)*2

    # marginalize over only one parameter; we now have two free parameters
    delta_chi2_68 = scipy.stats.chi2.ppf(pct68, 2)
    delta_chi2_95 = scipy.stats.chi2.ppf(pct95, 2)
    delta_chi2_997 = scipy.stats.chi2.ppf(pct997, 2)

    ax = fig.add_subplot(4,4,plotinds[(ind1,ind2)])
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
    ax.set_xlabel("{0} Value".format(par1.parname), fontsize=fontsize)
    ax.set_ylabel("{0} Value".format(par2.parname), fontsize=fontsize)
    pl.setp(ax.get_xticklabels(), rotation=90, fontsize=fontsize)
    pl.setp(ax.get_yticklabels(), fontsize=fontsize)

    ax.set_xlim(xmin, xmax,)
    ax.set_ylim(ymin, ymax,)

fig.tight_layout()

fig.subplots_adjust(hspace=0.7, wspace=0.7)
fig.savefig("ammonia_LTE_restricteddelta_error_estimate_demonstration.pdf")
