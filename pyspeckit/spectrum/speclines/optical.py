"""
Storage for optical spectral line information.
"""
from __future__ import print_function

import numpy as np

def hydrogen(nu,nl, vacuum=True):
    """
    Compute the rest wavelength of Hydrogen recombination lines in angstroms
    """
    rydberg = 10973731.6 # m^-1
    protontoelectron = 1836.15266 # ratio

    lvac = 1.0/rydberg * 1./(1/float(nl)**2 - 1/float(nu)**2) * 1e10 * (1.0+1.0/protontoelectron)

    if not vacuum:
        import ref_index
        return ref_index.vac2air(lvac/10)*10
    else:
        return lvac

# Format: name, units, vacuum?, display name
lines = {
         "H_alpha": [6564.614, 'Angstrom', True, r'$\mathrm{H}\alpha$'],
         "H_beta":  [4862.721, 'Angstrom', True, r'$\mathrm{H}\beta$'],
         "OIIIa":   [4960.295, 'Angstrom', True, r'$[\mathrm{OIII}]\lambda 4959\AA$'],
         "OIIIb":   [5008.239, 'Angstrom', True, r'$[\mathrm{OIII}]\lambda 5007\AA$'],
         "NIIa":    [6549.860, 'Angstrom', True, r'$[\mathrm{NII}]\lambda 6549\AA$'],
         "NIIb":    [6585.270, 'Angstrom', True, r'$[\mathrm{NII}]\lambda 6585\AA$'],
         "SIIa":    [6718.290, 'Angstrom', True, r'$[\mathrm{SII}]\lambda 6718\AA$'],
         "SIIb":    [6732.680, 'Angstrom', True, r'$[\mathrm{SII}]\lambda 6732\AA$'],
         "OI":      [6300.304, 'Angstrom', True, r'$[\mathrm{OI}]\lambda 6300\AA$'],
         "OII":     [3727.319, 'Angstrom', True, r'$[\mathrm{OII}]\lambda 3727\AA$'],
         "NeIII":   [3868.760, 'Angstrom', True, r'$[\mathrm{OII}]\lambda 3869\AA$']
}

def get_optical_lines():
    for i in range(3, 7):
        name = 'H_%i-2' % i
        wavelength = hydrogen(i, 2)
        lines[name] = [wavelength, 'Angstrom', True, name]

    xarr = []
    for key in lines.keys(): 
        xarr.append(lines[key][0])
    xarr = np.array(xarr)

    indx = np.argsort(xarr)
    xarr = np.sort(xarr)

    name = []
    keys = list(lines.keys())
    for i, key in enumerate(keys): 
        name.append(keys[indx[i]])
    name = np.array(name)

    xunits = []
    xvac = []
    dname = []
    for i, nombre in enumerate(name): 
        xunits.append(lines[nombre][1])
        xvac.append(lines[nombre][2])
        dname.append(lines[nombre][3])

    xunits = np.array(xunits)
    xvac = np.array(xvac)
    dname = np.array(dname)

    optical_lines = {'name': name, 'xarr': xarr, 'xunits': xunits, 'xvac': xvac, 'dname': dname}

    return optical_lines
