"""
LTE Molecule Modeling Tool
==========================

Uses astroquery to obtain molecular parameters.

Equations are based on Mangum & Shirley 2015 (2015PASP..127..266M)

Module API
^^^^^^^^^^
"""
from __future__ import print_function
import numpy as np
from astropy import units as u
from astropy import constants
from astropy import log
from .model import SpectralModel

kb_cgs = constants.k_B.cgs.value
h_cgs = constants.h.cgs.value
eightpicubed = 8 * np.pi**3
threehc = 3 * constants.h.cgs * constants.c.cgs
hoverk_cgs = (h_cgs/kb_cgs)
c_cgs = constants.c.cgs.value
ckms = constants.c.to(u.km/u.s).value

def line_tau(tex, total_column, partition_function, degeneracy, frequency,
             energy_upper, einstein_A=None):
    """
    Given the excitation temperature of the state, total column density of the
    molecule, the partition function, the degeneracy of the state, the
    frequency of the state, and the upper-state energy level, return the optical
    depth of that transition.

    This is a helper function for the LTE molecule calculations.  It implements
    the equations

    .. math::

        \\int \\tau_\\nu d\\nu = \\frac{c^2}{8 \\pi \\nu^2} A_{ij} N_u
                    \\left[ \\exp\\left( \\frac{h \\nu}{k_B T_{ex}} \\right)  - 1 \\right]

    or

    .. math::

        \\tau_\\nu = \\frac{c^2}{8 \\pi \\nu^2} A_{ij} N_u \\phi_\\nu
                    \\left[ \\exp\\left( \\frac{h \\nu}{k_B T_{ex}} \\right)  - 1 \\right]

    where

    .. math::
        N_{u} = N_{tot} \\frac{g_u}{Q} \\exp\\left(\\frac{-E_u}{k_B T_{ex}} \\right)

    based on Equations 11 and 29 of Mangum & Shirley 2015 (2015PASP..127..266M)

    The return value is therefore

    .. math::

        \\tau_\\nu / \\phi_\\nu


    The line profile function is, for a Gaussian, given by eqn A1:

    .. math::

        \\phi_\\nu = \\frac{1}{\\sqrt{2 \\pi} \\sigma}
        \\exp \\left[ -\\frac{(\\nu-\\nu_0)^2}{2 \\sigma^2} \\right]

    """
    # don't use dipole moment, because there are extra hidden dependencies

    assert frequency.unit.is_equivalent(u.Hz)
    assert energy_upper.unit.is_equivalent(u.erg)
    assert total_column.unit.is_equivalent(u.cm**-2)
    assert tex.unit.is_equivalent(u.K)

    N_upper = (total_column * degeneracy / partition_function *
               np.exp(-energy_upper / (constants.k_B * tex)))

    # equation 29 in Mangum 2015
    #taudnu = (eightpicubed * frequency * dipole_moment**2 / threehc *
    #             (np.exp(frequency*h_cgs/(kb_cgs*tex))-1) * N_upper)
    assert einstein_A.unit.is_equivalent(u.Hz)
    taudnu = ((constants.c**2/(8*np.pi*frequency**2) * einstein_A * N_upper)*
              (np.exp(frequency*constants.h/(constants.k_B*tex))-1))

    return taudnu.decompose()

# Deprecated version of the above
# def line_tau_nonquantum(tex, total_column, partition_function, degeneracy,
#                         frequency, energy_upper, SijMu2=None, molwt=None):
#
#     assert frequency.unit.is_equivalent(u.Hz)
#     assert energy_upper.unit.is_equivalent(u.erg)
#     assert total_column.unit.is_equivalent(u.cm**-2)
#     assert tex.unit.is_equivalent(u.K)
#
#     energy_lower = energy_upper - frequency*constants.h
#     #N_lower = (total_column * degeneracy / partition_function *
#     #           np.exp(-energy_lower / (constants.k_B * tex)))
#
#     # http://www.cv.nrao.edu/php/splat/OSU_Splat.html
#     assert SijMu2.unit.is_equivalent(u.debye**2)
#     amu = u.Da
#     assert molwt.unit.is_equivalent(amu)
#     C1 = 54.5953 * u.nm**2 * u.K**0.5 / amu**0.5 / u.debye**2
#     C2 = 4.799237e-5 * u.K / u.MHz
#     C3 = (1.43877506 * u.K / ((1*u.cm).to(u.Hz, u.spectral()) * constants.h)).to(u.K/u.erg)
#     #C3 = (constants.h / constants.k_B).to(u.K/u.erg)
#     tau = total_column/partition_function * C1 * (molwt/tex)**0.5 * (1-np.exp(-C2*frequency/tex)) * SijMu2 * np.exp(-C3*energy_lower/tex)
#
#     return tau.decompose()

def line_tau_cgs(tex, total_column, partition_function, degeneracy, frequency,
                 energy_upper, einstein_A):
    """
    Given the excitation temperature of the state, total column density of the
    molecule, the partition function, the degeneracy of the state, the
    frequency of the state, and the upper-state energy level, return the optical
    depth of that transition.

    Unlike :func:`line_tau`, this function requires inputs in CGS units.

    This is a helper function for the LTE molecule calculations.  It implements
    the equations

    .. math::

        \\int \\tau_\\nu d\\nu = \\frac{c^2}{8 \\pi \\nu^2} A_{ij} N_u
                    \\left[ \\exp\\left( \\frac{h \\nu}{k_B T_{ex}} \\right)  - 1 \\right]

    or

    .. math::

        \\tau_\\nu = \\frac{c^2}{8 \\pi \\nu^2} A_{ij} N_u \\phi_\\nu
                    \\left[ \\exp\\left( \\frac{h \\nu}{k_B T_{ex}} \\right)  - 1 \\right]

    where

    .. math::
        N_{u} = N_{tot} \\frac{g_u}{Q} \\exp\\left(\\frac{-E_u}{k_B T_{ex}} \\right)

    based on Equations 11 and 29 of Mangum & Shirley 2015 (2015PASP..127..266M)

    The return value is therefore

    .. math::

        \\tau_\\nu / \\phi_\\nu


    The line profile function is, for a Gaussian, given by eqn A1:

    .. math::

        \\phi_\\nu = \\frac{1}{\\sqrt{2 \\pi} \\sigma}
        \\exp \\left[ -\\frac{(\\nu-\\nu_0)^2}{2 \\sigma^2} \\right]

    """

    N_upper = (total_column * degeneracy / partition_function *
               np.exp(-energy_upper / (kb_cgs * tex)))

    # equation 29 in Mangum 2015
    #taudnu = (eightpicubed * frequency * dipole_moment**2 / threehc *
    #             (np.exp(frequency*h_cgs/(kb_cgs*tex))-1) * N_upper)

    # substitute eqn 11:
    # Aij = 64 pi^4 nu^3 / (3 h c^3) |mu ul|^2
    # becomes
    # dipole_moment^2 = 3 h c^3 Aij / ( 64 pi^4 nu^3 )

    taudnu = ((c_cgs**2/(8*np.pi*frequency**2) * einstein_A * N_upper)*
              (np.exp(frequency*h_cgs/(kb_cgs*tex))-1))

    return taudnu

def Jnu(nu, T):
    """RJ equivalent temperature (MS15 eqn 24)"""
    return constants.h*nu/constants.k_B / (np.exp(constants.h*nu/(constants.k_B*T))-1)

def Jnu_cgs(nu, T):
    """RJ equivalent temperature (MS15 eqn 24)
    (use cgs constants for speed)
    """
    return hoverk_cgs*nu / (np.exp(hoverk_cgs*nu/T)-1)

def line_brightness(tex, dnu, frequency, tbg=2.73*u.K, *args, **kwargs):
    tau = line_tau(tex=tex, frequency=frequency, *args, **kwargs) / dnu
    tau = tau.decompose()
    assert tau.unit == u.dimensionless_unscaled
    return (Jnu(frequency, tex)-Jnu(frequency, tbg)).decompose() * (1 - np.exp(-tau))

def line_brightness_cgs(tex, dnu, frequency, tbg=2.73, *args, **kwargs):
    tau = line_tau(tex=tex, frequency=frequency, *args, **kwargs) / dnu
    return (Jnu(frequency, tex)-Jnu(frequency, tbg)) * (1 - np.exp(-tau))

def get_molecular_parameters(molecule_name, tex=50, fmin=1*u.GHz, fmax=1*u.THz,
                             catalog='JPL',**kwargs):
    """
    Get the molecular parameters for a molecule from the JPL or CDMS catalog

    (this version should, in principle, be entirely self-consistent)

    Parameters
    ----------
    molecule_name : string
        The string name of the molecule (normal name, like CH3OH or CH3CH2OH,
        but it has to match the JPL catalog spec)
    tex : float
        Optional excitation temperature (basically checks if the partition
        function calculator works)
    catalog : 'JPL' or 'CDMS'
        Which catalog to pull from
    fmin : quantity with frequency units
    fmax : quantity with frequency units
        The minimum and maximum frequency to search over

    Examples
    --------
    >>> from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters
    >>> freqs, aij, deg, EU, partfunc = get_molecular_parameters('CH2CHCN',
    ...                                                          fmin=220*u.GHz,
    ...                                                          fmax=222*u.GHz,
                                                                )
    >>> freqs, aij, deg, EU, partfunc = get_molecular_parameters('CH3OH',
    ...                                                          fmin=90*u.GHz,
    ...                                                          fmax=100*u.GHz)
    """
    if catalog == 'JPL':
        from astroquery.jplspec import JPLSpec as QueryTool
    elif catalog == 'CDMS':
        from astroquery.linelists.cdms import CDMS as QueryTool
    else:
        raise ValueError("Invalid catalog specification")

    speciestab = QueryTool.get_species_table()
    if 'NAME' in speciestab.colnames:
        jpltable = speciestab[speciestab['NAME'] == molecule_name]
    elif 'molecule' in speciestab.colnames:
        jpltable = speciestab[speciestab['molecule'] == molecule_name]
    else:
        raise ValueError(f"Did not find NAME or molecule in table columns: {speciestab.colnames}")
    if len(jpltable) != 1:
        raise ValueError(f"Too many or too few matches to {molecule_name}")

    jpltbl = QueryTool.query_lines(fmin, fmax, molecule=molecule_name,
                                   parse_name_locally=True)

    def partfunc(tem):
        """
        interpolate the partition function
        WARNING: this can be very wrong
        """
        tem = u.Quantity(tem, u.K).value
        tems = np.array(jpltable.meta['Temperature (K)'])
        keys = [k for k in jpltable.keys() if 'q' in k.lower()]
        logQs = jpltable[keys]
        logQs = np.array(list(logQs[0]))
        inds = np.argsort(tems)
        #logQ = np.interp(tem, tems[inds], logQs[inds])
        # linear interpolation is appropriate; Q is linear with T... for some cases...
        # it's a safer interpolation, anyway.
        # to get a better solution, you can fit a functional form as shown in the
        # JPLSpec docs, but that is... left as an exercise.
        # (we can test the degree of deviation there)
        linQ = np.interp(tem, tems[inds], 10**logQs[inds])
        return linQ

    freqs = jpltbl['FREQ'].quantity
    freq_MHz = freqs.to(u.MHz).value
    deg = np.array(jpltbl['GUP'])
    EL = jpltbl['ELO'].quantity.to(u.erg, u.spectral())
    dE = freqs.to(u.erg, u.spectral())
    EU = EL + dE

    # need elower, eupper in inverse centimeter units
    elower_icm = jpltbl['ELO'].quantity.to(u.cm**-1).value
    eupper_icm = elower_icm + (freqs.to(u.cm**-1, u.spectral()).value)

    # from Brett McGuire https://github.com/bmcguir2/simulate_lte/blob/1f3f7c666946bc88c8d83c53389556a4c75c2bbd/simulate_lte.py#L2580-L2587

    # LGINT: Base 10 logarithm of the integrated intensity in units of nm2 ·MHz at 300 K.
    # (See Section 3 for conversions to other units.)
    # see also https://cdms.astro.uni-koeln.de/classic/predictions/description.html#description
    CT = 300
    logint = np.array(jpltbl['LGINT']) # this should just be a number
    #from CDMS website
    sijmu = (np.exp(np.float64(-(elower_icm/0.695)/CT)) - np.exp(np.float64(-(eupper_icm/0.695)/CT)))**(-1) * ((10**logint)/freq_MHz) * (24025.120666) * partfunc(CT)

    #aij formula from CDMS.  Verfied it matched spalatalogue's values
    aij = 1.16395 * 10**(-20) * freq_MHz**3 * sijmu / deg

    # we want logA for consistency with use in generate_model below
    aij = np.log10(aij)
    EU = EU.to(u.erg).value

    ok = np.isfinite(aij) & np.isfinite(EU) & np.isfinite(deg) & np.isfinite(freqs)

    return freqs[ok], aij[ok], deg[ok], EU[ok], partfunc



def generate_model(xarr, vcen, width, tex, column,
                   freqs, aij, deg, EU, partfunc,
                   background=None, tbg=2.73,
                   get_tau=False,
                   get_tau_sticks=False
                  ):
    """
    Model Generator

    Parameters
    ----------
    xarr : Quantity array [Hz]
        The X-axis (frequency)
    vcen : Quantity [km/s]
        The central velocity
    width : Quantity [km/s]
        The 1-sigma Gaussian line width [not FWHM]
    tex : Quantity [K]
        Excitation temperature
    column : Quantity [cm^-2]
        The total column density of the molecule.
        If specified as a value < 30, it will be assumed to be the log10 of the
        column density.
    freqs : Quantity array [Hz]
        The central frequency of the lines to model.
        Subsequent parameters, aij, deg, EU, must have
        the same shape as `freqs`.
    aij : array [log s^-1]
        The Einstein A coefficients, assumed to be in units
        of log10 of the rate in 1/s
    deg : array
        The transition degeneracy
    EU : Quantity array [erg]
        The upper state energy in ergs
    partfunc : function
        The partition function.  Can also be specified as an array of values in
        Kelvin.
    background : None
        An additional background field to include in the radiative transfer
        model.  Must be in Kelvin.
        [Presently not correctly implemented]
    tbg : Quantity [K]
        The background brightness temperature, defaulting to 2.73 for the CMB
        temperature.  This background will be treated as uniform with frequency
        and subtracted.
    get_tau : bool, default False
        If specified, the optical depth in each frequency bin will be returned
    get_tau_sticks : bool, default False
        If specified, the optical depth in each transition will be returned

    Examples
    --------
    Example case to produce a model::

        from pyspeckit.spectrum.models.lte_molecule import get_molecular_parameters, generate_model, generate_fitter
        freqs, aij, deg, eu, partfunc = get_molecular_parameters('CH3OH')
        def modfunc(xarr, vcen, width, tex, column):
            return generate_model(xarr, vcen, width, tex, column, freqs=freqs, aij=aij,
                                  deg=deg, EU=EU, partfunc=partfunc)

        fitter = generate_fitter(modfunc, name="CH3OH")
    """

    if hasattr(tex,'unit'):
        tex = tex.value
    # to allow for multi-excitation-temperature, we need a multizone model
    # the radiative transfer equation below is hard-coded to be single-zone

    if hasattr(column, 'unit'):
        column = column.value
    if np.isscalar(column):
        column = np.ones(len(freqs), dtype='float') * column
    else:
        assert len(column) == len(freqs)
        column = column.copy() # we're doing inplace modification below
    # assume low numbers are meant to be exponents
    low_col = column < 30
    log.debug(f"Found {low_col.sum()} entries that needed exponentiation")
    if any(low_col):
        column[low_col] = 10**column[low_col]

    if hasattr(tbg,'unit'):
        tbg = tbg.to(u.K).value
    if hasattr(vcen, 'unit'):
        vcen = vcen.to(u.km/u.s).value
    if hasattr(width, 'unit'):
        width = width.to(u.km/u.s).value

    #velo = xarr.to(u.km/u.s, equiv).value
    freq = xarr.to(u.Hz).value # same unit as nu below
    model = np.zeros_like(xarr).value

    # splatalogue can report bad frequencies as zero
    OK = (freqs.value != 0) & np.isfinite(aij) & np.isfinite(EU)

    freqs_ = freqs.to(u.Hz).value

    if callable(partfunc):
        Qs = partfunc(tex)
        if np.isscalar(Qs):
            Qs = np.ones(len(freqs), dtype='float') * Qs
    else:
        Qs = partfunc
        assert len(Qs) == len(freqs)

    model_tau = np.zeros_like(freq)
    tau_sticks = []

    for logA, gg, restfreq, eu, nt, Q in zip(aij[OK], deg[OK], freqs_[OK],
                                             EU[OK], column[OK], Qs[OK]):
        tau_over_phi = line_tau_cgs(tex=tex, total_column=nt,
                                    partition_function=Q, degeneracy=gg,
                                    frequency=restfreq, energy_upper=eu,
                                    einstein_A=10**logA)
        # commented out for performance log.debug(f"A={logA}, g={gg}, nu={restfreq}, eu={eu}, col={nt}, Q={Q}, tau/phi={tau_over_phi}")
        width_dnu = width / ckms * restfreq

        phi_nu = (
            ((2*np.pi)**0.5 * width_dnu)**-1 *
            np.exp(-(freq-(1-vcen/ckms)*restfreq)**2/(2*width_dnu**2)))

        tau_profile = (tau_over_phi * phi_nu)

        model_tau += tau_profile
        tau_sticks.append(tau_over_phi)

    if get_tau:
        return model_tau
    if get_tau_sticks:
        return tau_sticks

    jnu_line = Jnu_cgs(freq, tex)
    jnu_bg = Jnu_cgs(freq, tbg)
    jnu = (jnu_line-jnu_bg)

    # this is the same as below, but laid out more explicitly.  This form of the
    # equation implicity subtracts of a uniform-with-frequency background
    # model = jnu_line*(1-np.exp(-model_tau)) + jnu_bg*(np.exp(-model_tau)) - jnu_bg
    model = jnu*(1-np.exp(-model_tau))

    if background is not None:
        raise NotImplementedError
        return background-model
    return model


def generate_fitter(model_func, name):
    """
    Generator for hnco fitter class
    """

    myclass = SpectralModel(model_func, 4,
            parnames=['shift','width','tex','column'],
            parlimited=[(False,False),(True,False),(True,False),(True,False)],
            parlimits=[(0,0), (0,0), (0,0),(0,0)],
            shortvarnames=(r'\Delta x',r'\sigma','T_{ex}','N'),
            centroid_par='shift',
            )
    myclass.__name__ = name

    return myclass


def nupper_of_kkms(kkms, freq, Aul, replace_bad=None):
    r"""
    Mangum & Shirley 2015 eqn 82 gives, for the optically thin, Rayleigh-Jeans,
    negligible background approximation:

    .. math::

        N_{tot} = (3 k) / (8 \pi^3 \nu S \mu^2 R_i)   (Q/g) \exp(E_u/k T_{ex}) \int(T_R/f dv)

    Eqn 31:

    .. math::

        N_{tot}/N_u = Q_{rot} / g_u \exp(E_u/k T_{ex})

        -> N_{tot} = N_u Q_{rot} / g_u \exp(E_u/k T_{ex})
        -> Nu = N_{tot} g / Q_{rot} \exp(-E_u / k T_{ex})

    To get Nu of an observed line, then:

    .. math::

        Nu Q_{rot} / g_u \exp(E_u/k T_{ex}) = (3 k) / (8 \pi^3 \nu S \mu^2 R_i)   (Q/g) \exp(E_u/k T_{ex}) \int(T_R/f dv)

    This term cancels:

    .. math::
        Q_{rot} / g_u \exp(E_u/k T_{ex})

    Leaving:

    .. math::

        N_u = (3 k) / (8 \pi^3 \nu S \mu^2 R_i)   \int(T_R/f dv)

    $\int(T_R/f dv)$ is the optically thin integrated intensity in K km/s
    dnu/nu = dv/c [doppler eqn], so to get $\int(T_R dnu)$, sub in $dv = c/\nu d\nu$

    .. math::


        N_u = (3 k c) / (8 \pi^3 \nu^2  S \mu^2 R_i)   \int(T_R/f d\nu)


    We then need to deal with the S \mu^2 R_i term.  We assume R_i = 1, since we
    are not measuring any hyperfine transitions (R_i is the hyperfine
    degeneracy; eqn 75)
    Equation 11:

    .. math::

        A_{ul} = 64 \pi^4 \nu^3 / (3 h c^3) |\mu_{ul}|^2

    Equation 62:

    .. math::

        |\mu_{ul}|^2 = S \mu^2

        -> S \mu^2 = (3 h c^3 A_{ul}) / (64 \pi^4 \nu^3)

    Plugging that in gives

    .. math::

        Nu = (3 k c) / (8 \pi^3 \nu^2  ((3 h c^3 A_{ul}) / (64 \pi^4 \nu^3)))   \int(T_R/f d\nu)
           = (3 k c 64 \pi^4 \nu^3) / (8 \pi^3 \nu^2 3 h c^3 A_{ul})            \int(T_R/f d\nu)
           = (8 \pi \nu k / (A_{ul} c^2 h)) \int(T_R/f d\nu)

    which is the equation implemented below.  We could also have left this in
    dv units by substituting du = nu/c dv:

    .. math::

           = (8 \pi \nu^2 k / (A_{ul} c^3 h)) \int (T_R/f dv)

    """

    if replace_bad:
        neg = kkms <= 0
        kkms[neg] = replace_bad

    freq = u.Quantity(freq, u.GHz)
    Aul = u.Quantity(Aul, u.Hz)
    kkms = u.Quantity(kkms, u.K*u.km/u.s)

    nline = 8 * np.pi * freq * constants.k_B / constants.h / Aul / constants.c**2

    # kelvin-hertz
    Khz = (kkms * (freq/constants.c)).to(u.K * u.MHz)

    return (nline * Khz).to(u.cm**-2)

def ntot_of_nupper(nupper, eupper, tex, Q_rot, degeneracy=1):
    """ Given an N_upper, E_upper, tex, Q_rot, and degeneracy for a single state, give N_tot

    Mangum & Shirley 2015 eqn 31

    .. math::

        N_{tot}/N_u = Q_{rot} / g_u \\exp(E_u/k T_{ex})

    Example:
        >>> import astropy.units as u
        >>> tex = 50*u.K
        >>> kkms = 100*u.K*u.km/u.s
        >>> from pyspeckit.spectrum.models import lte_molecule
        >>> freqs, aij, deg, EU, partfunc = lte_molecule.get_molecular_parameters(molecule_name='HNCO v=0', molecule_name_jpl='HNCO', fmin=87*u.GHz, fmax=88*u.GHz)
        >>> nupper = lte_molecule.nupper_of_kkms(kkms, freqs, 10**aij)
        >>> ntot = lte_molecule.ntot_of_nupper(nupper, EU*u.erg, tex, Q_rot=partfunc(tex), degeneracy=deg)
    """

    Ntot = nupper * (Q_rot/degeneracy) * np.exp(eupper / (constants.k_B*tex))

    return Ntot

# url = 'http://cdms.ph1.uni-koeln.de/cdms/tap/'
# rslt = requests.post(url+"/sync", data={'REQUEST':"doQuery", 'LANG': 'VSS2', 'FORMAT':'XSAMS', 'QUERY':"SELECT SPECIES WHERE MoleculeStoichiometricFormula='CH2O'"})

if __name__ == "__main__":
    # example

    # 303
    J = 3
    gI = 0.25
    gJ = 2*J+1
    gK = 1

    # 321 has same parameters for g

    ph2co = {'tex':18.75*u.K,
             'total_column': 1e12*u.cm**-2,
             'partition_function': 44.6812, # splatalogue's 18.75
             'degeneracy': gI*gJ*gK,
             #'dipole_moment': 2.331e-18*u.esu*u.cm, #2.331*u.debye,
            }

    ph2co_303 = {
             'frequency': 218.22219*u.GHz,
             'energy_upper': kb_cgs*20.95582*u.K,
             'einstein_A': 10**-3.55007/u.s,
    }
    ph2co_303.update(ph2co)
    ph2co_303['dnu'] = (1*u.km/u.s/constants.c * ph2co_303['frequency'])

    ph2co_321 = {
             'frequency': 218.76007*u.GHz,
             'energy_upper': kb_cgs*68.11081*u.K,
             'einstein_A': 10**-3.80235/u.s,
    }
    ph2co_321.update(ph2co)
    ph2co_321['dnu'] = (1*u.km/u.s/constants.c * ph2co_321['frequency'])

    ph2co_322 = {
             'frequency': 218.47563*u.GHz,
             'energy_upper': kb_cgs*68.0937*u.K,
             'einstein_A': 10**-3.80373/u.s,
    }
    ph2co_322.update(ph2co)
    ph2co_322['dnu'] = (1*u.km/u.s/constants.c * ph2co_322['frequency'])

    print(("tau303 = {0}".format(line_tau(**ph2co_303))))
    print(("tau321 = {0}".format(line_tau(**ph2co_321))))
    print(("tau322 = {0}".format(line_tau(**ph2co_322))))
    print(("r303/r321 = {0}".format(line_brightness(**ph2co_321)/line_brightness(**ph2co_303))))
    print(("r303/r322 = {0}".format(line_brightness(**ph2co_322)/line_brightness(**ph2co_303))))

    # CDMS Q
    import requests
    import bs4
    url = 'http://cdms.ph1.uni-koeln.de/cdms/tap/'
    rslt = requests.post(url+"/sync", data={'REQUEST':"doQuery", 'LANG': 'VSS2', 'FORMAT':'XSAMS', 'QUERY':"SELECT SPECIES WHERE MoleculeStoichiometricFormula='CH2O'"})
    bb = bs4.BeautifulSoup(rslt.content, 'html5lib')
    h = [x for x in bb.findAll('molecule') if x.ordinarystructuralformula.value.text=='H2CO'][0]
    tem_, Q_ = h.partitionfunction.findAll('datalist')
    tem = [float(x) for x in tem_.text.split()]
    Q = [float(x) for x in Q_.text.split()]

    del ph2co_303['tex']
    del ph2co_303['partition_function']
    T_303 = np.array([line_brightness(tex=tex*u.K, partition_function=pf,
                                      **ph2co_303).value for tex,pf in
                      zip(tem,Q)])

    del ph2co_321['tex']
    del ph2co_321['partition_function']
    T_321 = np.array([line_brightness(tex=tex*u.K, partition_function=pf,
                                      **ph2co_321).value for tex,pf in
                      zip(tem,Q)])

    del ph2co_322['tex']
    del ph2co_322['partition_function']
    T_322 = np.array([line_brightness(tex=tex*u.K, partition_function=pf,
                                      **ph2co_322).value for tex,pf in
                      zip(tem,Q)])

    import pylab as pl

    pl.clf()
    pl.subplot(2,1,1)
    pl.plot(tem, T_321, label='$3_{2,1}-2_{2,0}$')
    pl.plot(tem, T_322, label='$3_{2,2}-2_{2,1}$')
    pl.plot(tem, T_303, label='$3_{0,3}-2_{0,2}$')
    pl.xlim(0,200)
    pl.subplot(2,1,2)
    pl.plot(tem, T_321/T_303, label='321/303')
    pl.plot(tem, T_322/T_303, label='322/303')
    pl.xlim(0,200)

    pl.draw(); pl.show()
