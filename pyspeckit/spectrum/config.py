"""
config.py

Author: Jordan Mirocha
Affiliation: University of Colorado at Boulder
Created on 2011-03-08.

Description: Modeled after config.py from yt.  Returns dictionary spcfg.cfg.

Notes: Will look in ~/.pyspeckit for file named config.  There is an example config file in the examples directory.
  
"""

import os, inspect
from functools import wraps

cfgDefaults = dict(
    color = 'k',
    composite_fit_color = 'r',
    component_fit_color = 'blue',
    lw = 0.5,
    composite_lw = 0.75,
    component_lw = 0.75,
    show_components = False,
    annotate = True,
    interactive = False,
    autorefresh = False,
    silent = True,
    debug = False
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
                
                if thisline[2] in ['True', 1]: return_dict[thisline[0]] = True
                elif thisline[2] in ['False', 0]: return_dict[thisline[0]] = False
                elif thisline[2] == 'None': return_dict[thisline[0]] = None
                elif thisline[2].isalpha(): return_dict[thisline[0]] = str(thisline[2])
                else: return_dict[thisline[0]] = float(thisline[2])
    	            	        
    	    self.cfg = return_dict
    	else: self.cfg = cfgDefaults
            	    	
__fn = os.path.expanduser("~/.pyspeckit/config")
if os.path.exists(__fn): 
    mycfg = ConfigParser(__fn).cfg
else:
    mycfg = ConfigParser().cfg
   
def ConfigDescriptor(f):    
                                
    def decorator(self, *args, **kwargs):  
        
        # inspect.getargspec will tell us the names of all arguments and their default values
        # later we'll have to be more careful - all_defs only makes entries for arguments that actually
        # have default values  
        all_args, all_vars, all_keys, all_defs = inspect.getargspec(f)                
        all_args.pop(0) # pop self
        
        # Construct dictionary containing all of f's keyword arguments                                   
        argsdefs = {}
        for i, arg in enumerate(all_args):
            argsdefs[arg] = all_defs[i]
                                                                   
        # Include these in our new_kwargs dictionary                                                                                                     
        new_kwargs = argsdefs
        
        # Read in config file and replace keyword arguments that have been defined in it
        for arg in new_kwargs:
            if mycfg.has_key(arg): new_kwargs[arg] = mycfg[arg]
                
        # If we've changed anything on call, reflect this in new_kwargs
        for arg in kwargs:
            #try: 
            #    if kwargs[arg] != argsdefs[arg]: new_kwargs[arg] = kwargs[arg]
            #except KeyError: 
            new_kwargs[arg] = kwargs[arg]       
                                                                                                                                       
        f(self, *args, **new_kwargs)
            
    return decorator      

