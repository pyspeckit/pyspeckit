"""
config.py

Author: Jordan Mirocha
Affiliation: University of Colorado at Boulder
Created on 2011-03-08.

Description: Modeled after config.py from yt.  Returns dictionary spcfg.cfg.

Notes: Will look in ~/.pyspeckit for file named config.

To do: Replace contents of confit with python dictionary?
  
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
    show_components = 0,
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
        all_args, all_vars, all_keys, all_defs = inspect.getargspec(f)                
        all_args.pop(0) # pop self
                                          
        argsdefs = {}
        for i, arg in enumerate(all_args):
            argsdefs[arg] = all_defs[i]
                                                                   
        # Include all supplied keyword arguments                                                                                                               
        new_kwargs = argsdefs
        
        # Read in config file
        for arg in new_kwargs:
            if mycfg.has_key(arg): new_kwargs[arg] = mycfg[arg]
                
        # Has anything changed?
        if f.__name__ == '__call__':
            for arg in kwargs:
                try: 
                    if kwargs[arg] != argsdefs[arg]: new_kwargs[arg] = kwargs[arg]
                except KeyError: 
                    new_kwargs[arg] = kwargs[arg]
                
                       
                        
        # Add all those from config file
       # for i, default in enumerate(all_defs):
       #     if new_kwargs.has_key(all_args[i]): continue
       #     
       #     if mycfg.has_key(all_args[i]): new_kwargs[all_args[i]] = mycfg[all_args[i]]
       #     else: new_kwargs[all_args[i]] = default
            
            # Only update parameter value if it isn't equal to its default...not sure this is absolutely correct
            #if kwargs.has_key(all_args[i]): 
            #    if kwargs[all_args[i]] != all_defs[i]:
            #        new_kwargs[all_args[i]] = kwargs[all_args[i]]          
                                                                                                               
        f(self, *args, **new_kwargs)
            
    return decorator      

