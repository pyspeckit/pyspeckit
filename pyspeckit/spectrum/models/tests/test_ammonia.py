import pytest
import numpy as np
from astropy import units as u
from ...classes import Spectrum, Spectra, units
from .. import ammonia, ammonia_constants

gg = [5., 2.8, 13.02918494, 0.0855, 14.85, 0.5]


def test_ammonia_parlimits():
    np.random.seed(0)

    sp = Spectrum(xarr=np.linspace(23.68,23.70)*u.GHz,
                  data=np.random.randn(50))

    sp.specfit(fittype='ammonia',
               guesses=gg,
              )

def test_ammonia_parlimits_fails():
    np.random.seed(0)

    sp = Spectrum(xarr=np.linspace(23.68,23.70)*u.GHz,
                  data=np.random.randn(50))

    with pytest.raises(ValueError) as ex:

        sp.specfit(fittype='ammonia',
                   guesses=gg,
                   limitedmin=[True, False, False, True, False, True],
                  )

    assert 'no such limit is set' in str(ex.value)

def make_synthspec(velo=np.linspace(-25, 25, 251), tkin=25, lte=True,
                   column=13, width=0.5, vcen=0.0, trot=None, tex=None,
                   lines=('oneone', 'twotwo', 'fourfour'),
                   ):

    if lte:
        trot = tex = tkin

    velo = u.Quantity(velo, u.km/u.s)

    spectra_list = []
    for line in lines:
        cfrq = ammonia_constants.freq_dict[line]*u.Hz
        frq = velo.to(u.GHz, u.doppler_radio(cfrq))
        xarr = units.SpectroscopicAxis(frq, refX=cfrq)
        data = ammonia.ammonia(xarr, tex=tex, trot=trot, width=width,
                               ntot=column, xoff_v=vcen)
        sp = Spectrum(xarr=xarr, data=data)
        spectra_list.append(sp)

    return Spectra(spectra_list)

def make_synthspec_cold(velo=np.linspace(-25, 25, 251), tkin=25,
                        column=13, width=0.5, vcen=0.0,
                        lines=('oneone', 'twotwo', 'fourfour'),
                        ):

    velo = u.Quantity(velo, u.km/u.s)

    spectra_list = []
    for line in lines:
        cfrq = ammonia_constants.freq_dict[line]*u.Hz
        frq = velo.to(u.GHz, u.doppler_radio(cfrq))
        xarr = units.SpectroscopicAxis(frq, refX=cfrq)
        data = ammonia.cold_ammonia(xarr, tkin=tkin, width=width, ntot=column,
                                    xoff_v=vcen)
        sp = Spectrum(xarr=xarr, data=data)
        spectra_list.append(sp)

    return Spectra(spectra_list)

def test_self_fit():
    """
    Self-consistency check: make sure inputs are recovered
    """
    spc = make_synthspec()
    spc.specfit(fittype='ammonia', guesses=[23, 22, 13.1, 1, 0.5, 0],
                fixed=[False,False,False,False,False,True])

    np.testing.assert_almost_equal(spc.specfit.parinfo[0].value, 25, 6)
    np.testing.assert_almost_equal(spc.specfit.parinfo[1].value, 25, 3)
    np.testing.assert_almost_equal(spc.specfit.parinfo[2].value, 13, 1)
    np.testing.assert_almost_equal(spc.specfit.parinfo[3].value, 0.5, 7)
    np.testing.assert_almost_equal(spc.specfit.parinfo[4].value, 0.0, 3)

def test_cold_ammonia():

    # from RADEX:
    # R = Radex(species='p-nh3', column=1e13, collider_densities={'pH2':1e4}, temperature=20)
    # R(collider_densities={'ph2': 1e4}, temperature=25, column=1e13)[:12]
    # trot = (u.Quantity(tbl['upperstateenergy'][8]-tbl['upperstateenergy'][9], u.K) *
    #    np.log((tbl['upperlevelpop'][9] * R.upperlevel_statisticalweight[8]) /
    #           (tbl['upperlevelpop'][8] * R.upperlevel_statisticalweight[9]))**-1
    #    )
    #    = 22.80
    spc = make_synthspec(lte=False, tkin=None, tex=6.66, trot=17.776914063182385,
                         lines=['oneone','twotwo'])
    spc.specfit.Registry.add_fitter('cold_ammonia',ammonia.cold_ammonia_model(),6)
    spc.specfit(fittype='cold_ammonia', guesses=[23, 5, 13.1, 1, 0.5, 0],
                fixed=[False,False,False,False,False,True]
               )

    np.testing.assert_almost_equal(spc.specfit.parinfo['tkin0'].value, 19.8365, 4)
    np.testing.assert_almost_equal(spc.specfit.parinfo['tex0'].value, 6.66, 3)
    np.testing.assert_almost_equal(spc.specfit.parinfo['ntot0'].value, 13, 2)
    np.testing.assert_almost_equal(spc.specfit.parinfo['width0'].value, 0.5, 6)
    np.testing.assert_almost_equal(spc.specfit.parinfo['xoff_v0'].value, 0.0, 3)

    return spc

def test_cold_ammonia2():

    tex = {'oneone':6.789360524825584,
           'twotwo':6.533403488783305,
          }
    spc = make_synthspec(lte=False, tkin=None, tex=tex, trot=17.776914063182385,
                         lines=['oneone','twotwo'])
    spc.specfit.Registry.add_fitter('cold_ammonia',ammonia.cold_ammonia_model(),6)
    spc.specfit(fittype='cold_ammonia', guesses=[23, 5, 13.1, 1, 0.5, 0],
                fixed=[False,False,False,False,False,True]
               )

    np.testing.assert_almost_equal(spc.specfit.parinfo['tkin0'].value, 19.55, 2)
    np.testing.assert_almost_equal(spc.specfit.parinfo['tex0'].value, 6.79, 2)
    np.testing.assert_almost_equal(spc.specfit.parinfo['ntot0'].value, 13, 2)
    np.testing.assert_almost_equal(spc.specfit.parinfo['width0'].value, 0.5, 4)
    np.testing.assert_almost_equal(spc.specfit.parinfo['xoff_v0'].value, 0.0, 3)

    return spc

def test_ammonia_x2():
    kwargs = dict(lte=False, tkin=None, lines=('oneone', 'twotwo'))
    redhot = dict(trot=25, tex=10, column=13, width=1, vcen=1)
    bluecold = dict(trot=12, tex=7, column=13, width=.2, vcen=-1)
    redhot.update(kwargs)
    bluecold.update(kwargs)
    sp_redhot = make_synthspec(**redhot)
    sp_bluecold = make_synthspec(**bluecold)

    sp = sp_redhot.copy()
    sp.Registry.add_fitter('ammonia', ammonia.ammonia_model(), 6)
    sp.specfit.fitter = sp.specfit.Registry.multifitters['ammonia']
    sp.specfit.fitter.npeaks = 2
    sp.data = sp_redhot.data + sp_bluecold.data
    assert np.allclose(sp.data.max(), 0.301252, rtol=1e-5)

    unpack_pars = lambda d: [d[key] for key in ['trot', 'tex',
                             'column', 'width', 'vcen']]+[0]

    model = sp.specfit.get_full_model(pars=unpack_pars(redhot)
                                      +unpack_pars(bluecold))
    assert np.allclose(sp.data, model, rtol=1e-5)

def test_regression_179():
    # making a dummy spectrum
    xarr = units.SpectroscopicAxis(range(20), 'km/s')
    xarr.velocity_convention='radio'
    xarr.refX = ammonia_constants.freq_dict['oneone']*u.Hz

    # one zero-amplitude component and one with Tmb>0
    blank = [10., 2.7315, 14., 0.2, 0.0, 0.5]
    nonzero = [10., 5, 14., 0.2, 10.0, 0.5]

    c_nh3 = ammonia.cold_ammonia_model()

    # so far so good, the modeled amplitudes are okay
    assert c_nh3.n_modelfunc(blank)(xarr).max() == 0
    assert c_nh3.n_modelfunc(nonzero)(xarr).max() > 0

    # but I was expecting the two below to be the same
    assert c_nh3.n_modelfunc(blank+nonzero)(xarr).max() > 0
    assert c_nh3.n_modelfunc(nonzero+blank)(xarr).max() > 0

    sp = Spectrum(data=np.zeros_like(xarr.value), xarr=xarr)
    sp.Registry.add_fitter('cold_ammonia', ammonia.cold_ammonia_model(),6)
    sp.specfit.fitter = sp.specfit.Registry.multifitters['cold_ammonia']
    sp.specfit.fittype = 'cold_ammonia'
    assert sp.specfit.get_full_model(pars=nonzero+blank).max() > 0


"""
# debug work to make sure synthspectra look the same
from pyspeckit.spectrum.models.tests import test_ammonia
spc1 = test_ammonia.test_cold_ammonia()
spc2 = test_ammonia.test_cold_ammonia2()
spc1.plotter()
spc2.plotter(axis=spc1.plotter.axis, color='b', clear=False)

from pyspeckit.spectrum.models.tests import test_ammonia
from astropy import log
log.setLevel('DEBUG')
spc3 = test_ammonia.make_synthspec(lte=False, tkin=None, trot=17.8, tex=6.66, lines=['oneone','twotwo'])
spc4 = test_ammonia.make_synthspec(lte=False, tkin=None, trot=17.8, tex={'oneone':6.66, 'twotwo':6.66}, lines=['oneone','twotwo'])
spc5 = test_ammonia.make_synthspec(lte=False, tkin=None, trot=17.8, tex={'oneone':6.66, 'twotwo':6.66, 'fourfour':6.66}, lines=['oneone','twotwo'])
spc3.plotter(linewidth=2)
spc5.plotter(axis=spc3.plotter.axis, color='r', clear=False)
spc4.plotter(axis=spc3.plotter.axis, color='b', clear=False)
"""

def test_ammonia_guessing():

    mod = ammonia.ammonia_model()

    rslt = mod.parse_3par_guesses([1.0, 'cen', 'wid'])

    assert rslt == [3.73*2, 3.73, 15, 'wid', 'cen', 0.5]

    rslt = mod.parse_3par_guesses([1.0, 'cen1', 'wid1',
                                   2.0, 'cen2', 'wid2',
                                  ])

    assert rslt == [3.73*2, 3.73, 15, 'wid1', 'cen1', 0.5,
                    4.73*2, 4.73, 15, 'wid2', 'cen2', 0.5,
                   ]
