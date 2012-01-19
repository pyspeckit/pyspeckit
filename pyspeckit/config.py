"""
==========================
PySpecKit config system
==========================

Create decorator used to modify inputs to various __call__ methods
to reflect preferences we set in ~/.pyspeckit/config.

To see what values exist in config file, do:
    from config import mycfg
    mycfg.keys()
    
To decorate a new __call__ method, do:
    from config import ConfigDescriptor as cfgdec
    
    @cfgdec
    def __call__(self, **kwargs):
        pass    # do something!

.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>
.. moduleauthor:: Jordan Mirocha <mirochaj@gmail.com>
"""

import os, inspect

class dotdictify(dict):
    """
    Class grabbed from
    http://stackoverflow.com/questions/3031219/python-recursively-access-dict-via-attributes-as-well-as-index-access
    to allow "dot" access to dictionary keys
    """
    marker = object()
    def __init__(self, value=None):
        if value is None:
            pass
        elif isinstance(value, dict):
            for key in value:
                self.__setitem__(key, value[key])
        else:
            raise TypeError, 'expected dict'

    def __setitem__(self, key, value):
        if isinstance(value, dict) and not isinstance(value, dotdictify):
            value = dotdictify(value)
        dict.__setitem__(self, key, value)

    def __getitem__(self, key):
        found = self.get(key, dotdictify.marker)
        if found is dotdictify.marker:
            found = dotdictify()
            dict.__setitem__(self, key, found)
        return found

    __setattr__ = __setitem__
    __getattr__ = __getitem__


cfgDefaults = dotdictify(dict(
    color = 'k',
    composite_fit_color = 'r',
    component_fit_color = 'blue',
    baseline_fit_color = 'orange',
    lw = 0.5,
    composite_lw = 0.75,
    component_lw = 0.75,
    show_components = False,
    autoannotate = True,
    interactive = False,
    autorefresh = True,
    silent = True,
    debug = False,
    WARN = True,
    ))

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
    	            	        
    	    self.cfg = dotdictify(return_dict)
    	else: 
    	    self.cfg = dotdictify(cfgDefaults)
            	    	
__fn = os.path.expanduser("~/.pyspeckit/config")
if os.path.exists(__fn): 
    mycfg = dotdictify(ConfigParser(__fn).cfg)
else:
    mycfg = dotdictify(ConfigParser().cfg)
   
def ConfigDescriptor(f):    
                                
    def decorator(self, *args, **kwargs):
        """
        This is our decorator function, used to modify the inputs of __call__
        methods to reflect preferences we set in config file.
                
        Notes:
        inspect.getargspec will tell us the names of all arguments and their default values.
        Later we'll have to be more careful - all_defs only makes entries for arguments that actually
        have default values          
        """  

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
            if arg in mycfg: new_kwargs[arg] = mycfg[arg]
                
        # If we've changed anything on call, reflect this in new_kwargs
        for arg in kwargs:
            new_kwargs[arg] = kwargs[arg]       
                                                                                                                                       
        f(self, *args, **new_kwargs)
            
    # documentation should be passed on, else sphinx doesn't work and the user can't access the docs
    decorator.__doc__ = f.__doc__
    decorator.__defaults__ = f.__defaults__
    decorator.__repr__ = f.__repr__
    return decorator      

