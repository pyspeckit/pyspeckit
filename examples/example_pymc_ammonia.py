import numpy as np
from astropy import units as u
import itertools
from operator import itemgetter
import pyspeckit
import scipy.stats

from pyspeckit.spectrum.models import ammonia, ammonia_constants

import pylab as pl
pl.close('all')
pl.figure(1).clf()

oneonefreq = ammonia_constants.freq_dict['oneone']
twotwofreq = ammonia_constants.freq_dict['twotwo']
# create an axis that covers the 1-1 and 2-2 inversion lines
xaxis = pyspeckit.units.SpectroscopicAxis(np.linspace(oneonefreq*(1-50/3e5),
                                                      twotwofreq*(1+50/3e5),
                                                      1000.),
                                          unit=u.Hz)
sigma = 2.
center = 0.
trot = 15.
ntot = 14.7
# Adopting an NLTE model: T_ex < T_rot (this case is better for the default fitter)
# pyradex.Radex(species='p-nh3', column=5e14, collider_densities={'h2':5000}, temperature=15)()[8:10]
tex = {'oneone':8.5, 'twotwo':5.5}
synth_data = ammonia.ammonia(xaxis, trot=trot, tex=tex, width=sigma,
                             xoff_v=center, ntot=ntot)

# Add noise
stddev = 0.1
noise = np.random.randn(xaxis.size)*stddev
error = stddev*np.ones_like(synth_data)
data = u.Quantity(noise+synth_data, u.K)

# create the spectrum object
sp = pyspeckit.Spectrum(data=data,
                        xarr=xaxis,
                        error=np.ones_like(synth_data)*stddev)
# fit it
# Setting limits is important for running the MCMC chains below (they set the
# priors)
sp.specfit(fittype='ammonia',
           # trot, tex, ntot, width, center, fortho
           guesses=[35,10,15,5,5,0],
           fixed=[False,False,False,False,False,True],
           minpars=[3,3,10,0,-10,0],
           maxpars=[80,80,20,10,10,1],
           limitedmin=[True]*6,
           limitedmax=[True]*6,
          )
# then get the pymc values
MCuniformpriors = sp.specfit.get_pymc()
MCuniformpriors.sample(10100,burn=100,tune_interval=250)

# MC vs least squares:
print(sp.specfit.parinfo)

print(MCuniformpriors.stats()['tex0'],
      MCuniformpriors.stats()['ntot0'])



# optional plotting
from mpl_plot_templates import pymc_plotting
import pylab
pylab.figure(1).clf()
pymc_plotting.hist2d(MCuniformpriors, 'tex0', 'ntot0',
                     bins=[25,25])
pylab.plot([tex['oneone']],[ntot],'k+',markersize=25)
pylab.axis([6,10.5,14.6,14.9])
pylab.savefig("tex_vs_ntot_pymc_example.pdf")


# Now do the same with emcee
emcee_ensemble = sp.specfit.get_emcee()
p0 = emcee_ensemble.p0 * (np.random.randn(*emcee_ensemble.p0.shape) / 50. + 1.0)
pos,logprob,state = emcee_ensemble.run_mcmc(p0,10000)

plotdict = {'tex0':emcee_ensemble.chain[:,1000:,1].ravel(),
            'ntot0':emcee_ensemble.chain[:,1000:,2].ravel()}
# one of the samplers went off the hook and got locked at high ntot
OK = (plotdict['ntot0'] < 20)
plotdict['tex0'] = plotdict['tex0'][OK]
plotdict['ntot0'] = plotdict['ntot0'][OK]

pylab.figure(2).clf()
pymc_plotting.hist2d(plotdict, 'tex0', 'ntot0', fignum=2, bins=[25,25],
                     clear=True)
pylab.plot([tex['oneone']],[ntot],'k+',markersize=25)
pylab.axis([6,10.5,14.6,14.9])
pylab.savefig("tex_vs_ntot_emcee_example.pdf")
