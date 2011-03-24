"""
config.py

Author: Jordan Mirocha
Affiliation: University of Colorado at Boulder
Created on 2011-03-08.

Description: Modeled after config.py from yt.  Returns dictionary spcfg.cfg.

Notes: Will look in ~/.pyspeckit for file named config.

Sample ~/.spectools/config file:

plot_color = k
fit_color = green
comp_color = blue
plot_lw = 0.5
fit_lw = 1.0
comp_lw = 1.0

To do: choose reasonable naming convention for these params.
  
"""

import os

cfgDefaults = dict(
    plot_color = 'k',
    fit_color = 'r',
    comp_color = 'blue',
    plot_lw = 0.5,
    fit_lw = 0.75,
    comp_lw = 0.75,
    show_components = 0,
    # FITS should be the default; it is the standard for most spectrographs
    data_format = 'fits', 
    annotate = True,
    logfile='SpectrumLogger.log',
    )

class ConfigParser:
    def __init__(self, fn = None):
        """
        Initialize cfg dictionary.
        """
        if fn is not None:
            f = open(fn)
            return_dict = cfgDefaults
            for line in f:
                if not line.strip(): continue
                thisline = line.split()
            
                if thisline[0][0] == '#': continue
                
                try: return_dict[thisline[0]] = float(thisline[2])
                except ValueError: return_dict[thisline[0]] = str(thisline[2])
    	        
    	    self.cfg = return_dict
    	else: self.cfg = cfgDefaults
    	
__fn = os.path.expanduser("~/.pyspeckit/config")
if os.path.exists(__fn): 
    spcfg = ConfigParser(__fn)
else:
    spcfg = ConfigParser()
