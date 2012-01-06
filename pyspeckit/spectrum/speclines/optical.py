"""
Storage for optical spectral line information.
"""

import numpy as np

# Format: name, units, vacuum?, display name
lines = {
         "H_alpha": [6564.614, 'Angstrom', True, r'$\mathrm{H}\alpha$'], \
         "H_beta":  [4862.721, 'Angstrom', True, r'$\mathrm{H}\beta$'], \
         "OIIIa":   [4960.295, 'Angstrom', True, r'$[\mathrm{OIII}]\lambda 4959\AA$'], \
         "OIIIb":   [5008.239, 'Angstrom', True, r'$[\mathrm{OIII}]\lambda 5007\AA$'], \
         "NIIa":    [6549.860, 'Angstrom', True, r'$[\mathrm{NII}]\lambda 6549\AA$'], \
         "NIIb":    [6585.270, 'Angstrom', True, r'$[\mathrm{NII}]\lambda 6585\AA$'], \
         "SIIa":    [6718.290, 'Angstrom', True, r'$[\mathrm{SII}]\lambda 6718\AA$'], \
         "SIIb":    [6732.680, 'Angstrom', True, r'$[\mathrm{SII}]\lambda 6732\AA$'], \
         "OI":      [6300.304, 'Angstrom', True, r'$[\mathrm{OI}]\lambda 6300\AA$'], \
         "OII":     [3727.319, 'Angstrom', True, r'$[\mathrm{OII}]\lambda 3727\AA$'], \
         "NeIII":   [3868.760, 'Angstrom', True, r'$[\mathrm{OII}]\lambda 3869\AA$']
}

xarr = []
for key in lines.keys(): 
    xarr.append(lines[key][0])
xarr = np.array(xarr)

indx = np.argsort(xarr)
xarr = np.sort(xarr)

name = []
for i, key in enumerate(lines.keys()): 
    name.append(lines.keys()[indx[i]])
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
