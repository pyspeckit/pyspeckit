import numpy as np
from astropy import units as u
from astropy import constants

kb_cgs = constants.k_B.cgs
h_cgs = constants.h.cgs
eightpicubed = 8 * np.pi**3
threehc = 3 * constants.h.cgs * constants.c.cgs


def line_tau(tex, total_column, partition_function, degeneracy, frequency,
             dipole_moment, energy_upper):

    N_upper = (total_column * degeneracy / partition_function *
               np.exp(-energy_upper / (kb_cgs * tex)))

    # equation 29 in Mangum 2015
    tauperdnu = (eightpicubed * frequency * dipole_moment**2 / threehc *
                 (np.exp(h_cgs*frequency/(kb_cgs*tex))-1) * N_upper)

    return tauperdnu.decompose()

def Jnu(nu, T):
    """RJ equivalent temperature (eqn 24)"""
    return constants.h*nu/constants.k_B / (np.exp(constants.h*nu/(constants.k_B*T))-1)

def line_brightness(tex, dnu, frequency, tbg=2.73*u.K, *args, **kwargs):
    tau = line_tau(tex=tex, frequency=frequency, *args, **kwargs) / dnu
    tau = tau.decompose()
    assert tau.unit == u.dimensionless_unscaled
    return (Jnu(frequency, tex)-Jnu(frequency, tbg)) * (1 - np.exp(-tau))

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
             'dipole_moment': 2.331e-18*u.esu*u.cm, #2.331*u.debye,
            }

    ph2co_303 = {
             'frequency': 218.22219*u.GHz,
             'energy_upper': kb_cgs*20.95582*u.K,
    }
    ph2co_303.update(ph2co)
    ph2co_303['dnu'] = (1*u.km/u.s/constants.c * ph2co_303['frequency'])

    ph2co_321 = {
             'frequency': 218.76007*u.GHz,
             'energy_upper': kb_cgs*68.11081*u.K,
    }
    ph2co_321.update(ph2co)
    ph2co_321['dnu'] = (1*u.km/u.s/constants.c * ph2co_321['frequency'])

    ph2co_322 = {
             'frequency': 218.47563*u.GHz,
             'energy_upper': kb_cgs*68.0937*u.K,
    }
    ph2co_322.update(ph2co)
    ph2co_322['dnu'] = (1*u.km/u.s/constants.c * ph2co_322['frequency'])

    print("tau303 = {0}".format(line_tau(**ph2co_303)))
    print("tau321 = {0}".format(line_tau(**ph2co_321)))
    print("tau322 = {0}".format(line_tau(**ph2co_322)))
    print("r303/r321 = {0}".format(line_brightness(**ph2co_321)/line_brightness(**ph2co_303)))
    print("r303/r322 = {0}".format(line_brightness(**ph2co_322)/line_brightness(**ph2co_303)))

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
