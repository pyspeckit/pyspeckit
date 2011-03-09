"""
speclines.py

Author: Jordan Mirocha
Affiliation: University of Colorado at Boulder
Created on 2011-03-09.

Description: Container for spectral line information.
"""

import numpy as np

# Common lines in the optical
opt_waves = np.array(
             [1033.30, 1215.67, 1239.42, 1305.53, 1335.52, 1399.8, 1545.86, 1640.4, 
              1665.85, 1857.4, 1908.27, 2326.0, 2439.5, 2800.32, 3346.79, 3426.85,
              3728.30, 3798.976, 3836.47, 3889.0, 3934.777, 3969.588, 4072.3, 
              4102.89, 4305.61, 4341.68, 4364.436, 4862.68, 4960.295, 5008.240,
              5176.7, 5895.6, 6302.046, 6365.536, 6549.86, 6564.61, 6585.27,
              6707.89, 6718.29, 6732.67])
opt_lines = np.array(
             ['OVI', 'Ly_alpha', 'NV', 'OI', 'CII', 'SiIV+OIV', 'CIV', 'HeII', 
              'OIII', 'AlIII', 'CIII', 'CII', 'NeIV', 'MgII', 'NeV', 'NeV', 'OII', 
              'H_theta', 'H_eta', 'HeI', 'K', 'H', 'SII', 'H_delta', 'G', 'H_gamma',
              'OIII', 'H_beta', 'OIII', 'OIII', 'Mg', 'Na', 'OI', 'OI', 'NII', 'H_alpha',
              'NII', 'Li', 'SII', 'SII'])
              
optical_lines = {'line': opt_lines, 'wavelength': opt_waves, 'units': 'Angstrom', 'vac': True}