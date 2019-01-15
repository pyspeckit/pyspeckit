import pyspeckit
import numpy as np
from pyspeckit.spectrum.models import inherited_voigtfitter

# technically, the voigt fitter works as a singlefitter (i.e., you can fit the
# background level and the peak simultaneously) 
# in practice, however, you need to fit the background independently except for
# gaussians.  I don't know why this is.

xarr = pyspeckit.spectrum.units.SpectroscopicAxis(np.linspace(-100,100,500),unit='km/s',refX=1e9,refX_unit='Hz')
VF = inherited_voigtfitter.voigt_fitter()
nvoigt = VF.n_modelfunc

synthdata = nvoigt([1,-30,6.5,0.5,0.5,35,1.5,6.5])(xarr) + np.random.randn(xarr.shape[0])/20.
sp1 = pyspeckit.Spectrum(xarr=xarr, data=synthdata,
                         error=np.ones(xarr.shape[0])/20., header={},)
sp1.plotter()
sp1.specfit(fittype='gaussian', guesses=[0.5,-25,3,0.2,40,5],
            composite_fit_color='b', clear=False, annotate=False)
sp1.specfit(fittype='lorentzian', guesses=[0.5,-25,3,0.2,40,5],
            composite_fit_color='g', clear=False, annotate=False)

sp1.specfit(fittype='voigt', guesses=[0.5,-30,2,2,0.5,45,2,2],
            tied=['','','','','','p[1]+65','',''],
            composite_fit_color='r',clear=False,annotate=True)
sp1.baseline(excludefit=True)
sp1.baseline.annotate()



# this approach doesn't work right now, but it will (there's a bug I'm working on)
# it's a lot more verbose, so it's kinda ugly, but it is (in principle) more flexible
parinfo = pyspeckit.parinfo.ParinfoList()
parinfo.append(pyspeckit.parinfo.Parinfo(parname='AMP',value=0.5))
parinfo.append(pyspeckit.parinfo.Parinfo(parname='VELO',value=-30))
parinfo.append(pyspeckit.parinfo.Parinfo(parname='GWIDTH',value=2))
parinfo.append(pyspeckit.parinfo.Parinfo(parname='LWIDTH',value=2))
parinfo.append(pyspeckit.parinfo.Parinfo(parname='AMP',value=0.5))
parinfo.append(pyspeckit.parinfo.Parinfo(parname='VELO',value=35, tied='p[1]+65'))
parinfo.append(pyspeckit.parinfo.Parinfo(parname='GWIDTH',value=2))
parinfo.append(pyspeckit.parinfo.Parinfo(parname='LWIDTH',value=2))

sp1.specfit(fittype='voigt', parinfo=parinfo, multifit=None,
            composite_fit_color='c', clear=False, annotate=True)

# Want to know the FWHM?
fwhm = sp1.specfit.fitter.analytic_fwhm()
real_fwhm = [inherited_voigtfitter.voigt_fwhm(6.5,0.5),
             inherited_voigtfitter.voigt_fwhm(1.5,6.5)]
print(fwhm,real_fwhm)
# Check that the real and analytic are close enough, given the noise
assert abs(fwhm[0]-real_fwhm[0]) < 1.5
