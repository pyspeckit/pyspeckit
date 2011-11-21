"""
==========================
PySpecKit HDF5 Reader
==========================

Routines for reading in spectra from HDF5 files.

Note: Current no routines for parsing HDF5 headers in classes.py.

.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>
.. moduleauthor:: Jordan Mirocha <mirochaj@gmail.com>
"""
try: 
    import h5py
    h5OK = True
except ImportError: 
    h5OK = False
    
import numpy as np
from .. import units

try: 
    import pyfits
    pyfitscheck=True
except ImportError: 
    pyfitscheck = False

def open_hdf5(filename, xaxkey = 'xarr', datakey = 'data', errkey = 'error'):
    """
    This reader expects three datasets to exist in the hdf5 file 'filename': 
    'xarr', 'data', and 'error', by default.  Can specify other dataset names.
    """
    
    if not h5OK: 
        print "WARNING: h5py not installed; cannot read hdf5 files."
    else:
        f = h5py.File(filename, 'r')
        
        try: xarr = f[xaxkey].value
        except KeyError: print 'Dataset \'%s\' not found.' % xaxkey
        
        try: data = f[datakey].value
        except KeyError: print 'Dataset \'%s\' not found.' % datakey
    
        try: error = f[errkey].value
        except KeyError: 
            print 'Dataset \'%s\' not found.  Assuming uniform errors.' % errkey
            error = np.ones_like(data)
            
        try: 
            xunits = f[xaxkey].attrs['units']
            if xunits == '': 
                xunits = 'unknown'
        except KeyError: 
            xunits = 'unknown'    
        
        try: 
            yunits = f[datakey].attrs['units']
            if yunits == '': 
                yunits = 'unknown'
        except KeyError: 
            yunits = 'unknown'    
            
        try: 
            xtype = f[xaxkey].attrs['type']
            if xtype == '': 
                xtype = xaxkey
        except KeyError: 
            xtype = xaxkey
        
        try: 
            ytype = f[datakey].attrs['type']
            if ytype == '': 
                ytype = datakey
        except KeyError: 
            ytype = datakey            
            
        XAxis = units.SpectroscopicAxis(xarr, xunits)    

        header = {'xunits': xunits, 'xtype': xtype, 'yunits': yunits, 'ytype': ytype}
        
        return data, error, XAxis, header
    
