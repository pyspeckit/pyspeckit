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
    #    np.log((tbl['upperlevelpop'][9] * R.statistical_weight[8]) /
    #           (tbl['upperlevelpop'][8] * R.statistical_weight[9]))**-1
    #    )
    #    = 22.80
    spc = make_synthspec(lte=False, tkin=None, tex=6.66, trot=17.776914063182385)
    spc.specfit.Registry.add_fitter('cold_ammonia',ammonia.cold_ammonia_model(),6)
    spc.specfit(fittype='cold_ammonia', guesses=[23, 5, 13.1, 1, 0.5, 0],
                fixed=[False,False,False,False,False,True]
               )

    np.testing.assert_almost_equal(spc.specfit.parinfo['tkin0'].value, 19.8365, 4)
    np.testing.assert_almost_equal(spc.specfit.parinfo['tex0'].value, 6.66, 3)
    np.testing.assert_almost_equal(spc.specfit.parinfo['ntot0'].value, 13, 2)
    np.testing.assert_almost_equal(spc.specfit.parinfo['width0'].value, 0.5, 7)
    np.testing.assert_almost_equal(spc.specfit.parinfo['xoff_v0'].value, 0.0, 4)

    return spc
