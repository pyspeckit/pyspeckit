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

def open_hdf5(filename):
    """
    This reader expects three datasets to exist in the hdf5 file 'filename': 'xarr', 'data', and 'error'.
    """
    
    if not h5OK: print "Cannot read hdf5 format - h5py import failed."
    else:
        f = h5py.File(filename, 'r')
        
        try: xarr = f['xarr'].value
        except KeyError: print 'Dataset \'xarr\' not found.'
        
        try: data = f['data'].value
        except KeyError: print 'Dataset \'data\' not found.'
    
        try: error = f['error'].value
        except KeyError: 
            print 'Dataset \'err\' not found.'
            error = np.ones_like(data)
            
        try: 
            xunits = f['xarr'].attrs['xunits']
            if xunits == '': xunits = 'unknown'
        except KeyError: xunits = 'unknown'    
            
        XAxis = units.SpectroscopicAxis(xarr,xunits)    

        if pyfitscheck:
            header = pyfits.Header()
            header.update('CUNIT1',xunits)
        else:
            header = {}
        
        return data,error,XAxis,header
    
