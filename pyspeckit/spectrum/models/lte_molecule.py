from __future__ import print_function
import numpy as np
from astropy import units as u
from astropy import constants

kb_cgs = constants.k_B.cgs.value
h_cgs = constants.h.cgs.value
eightpicubed = 8 * np.pi**3
threehc = 3 * constants.h.cgs * constants.c.cgs
hoverk_cgs = (h_cgs/kb_cgs)
c_cgs = constants.c.cgs.value


def line_tau(tex, total_column, partition_function, degeneracy, frequency,
             energy_upper, einstein_A=None):
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

def line_tau_nonquantum(tex, total_column, partition_function, degeneracy,
                        frequency, energy_upper, SijMu2=None, molwt=None):

    assert frequency.unit.is_equivalent(u.Hz)
    assert energy_upper.unit.is_equivalent(u.erg)
    assert total_column.unit.is_equivalent(u.cm**-2)
    assert tex.unit.is_equivalent(u.K)

    energy_lower = energy_upper - frequency*constants.h
    #N_lower = (total_column * degeneracy / partition_function *
    #           np.exp(-energy_lower / (constants.k_B * tex)))

    # http://www.cv.nrao.edu/php/splat/OSU_Splat.html
    assert SijMu2.unit.is_equivalent(u.debye**2)
    amu = u.Da
    assert molwt.unit.is_equivalent(amu)
    C1 = 54.5953 * u.nm**2 * u.K**0.5 / amu**0.5 / u.debye**2
    C2 = 4.799237e-5 * u.K / u.MHz
    C3 = (1.43877506 * u.K / ((1*u.cm).to(u.Hz, u.spectral()) * constants.h)).to(u.K/u.erg)
    #C3 = (constants.h / constants.k_B).to(u.K/u.erg)
    tau = total_column/partition_function * C1 * (molwt/tex)**0.5 * (1-np.exp(-C2*frequency/tex)) * SijMu2 * np.exp(-C3*energy_lower/tex)

    return tau.decompose()

def line_tau_cgs(tex, total_column, partition_function, degeneracy, frequency,
                 energy_upper, einstein_A):

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
    """RJ equivalent temperature (eqn 24)"""
    return constants.h*nu/constants.k_B / (np.exp(constants.h*nu/(constants.k_B*T))-1)

def Jnu_cgs(nu, T):
    """RJ equivalent temperature (eqn 24)
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

class lte_line_model_generator(object):
    def __init__(self, molecule_name, **kwargs):
        pass
        
    def line_profile(xarr,):
        pass

# requires vamdc branch of astroquery
def get_molecular_parameters(molecule_name, chem_re_flags=0):
    from astroquery.vamdc import load_species_table

    lut = load_species_table()
    species_id_dict = lut.find(molecule_name, flags=chem_re_flags)
    if len(species_id_dict) == 1:
        species_id = list(species_id_dict.values())[0]
    else:
        raise ValueError("Too many species matched: {0}"
                         .format(species_id_dict))
     
    request = r.Request(node=cdms)
    query_string = "SELECT ALL WHERE VAMDCSpeciesID='%s'" % species_id
    request.setquery(query_string)
    result = request.dorequest()
    Q = m.calculate_partitionfunction(result.data['States'],
                                      temperature=tex)[species_id]

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
