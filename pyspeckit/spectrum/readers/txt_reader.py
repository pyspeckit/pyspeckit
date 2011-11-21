"""
==========================
PySpecKit ASCII Reader
==========================

Routines for reading in ASCII format spectra.  If atPy is not installed,
will use a very simple routine for reading in the data.  

.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>
.. moduleauthor:: Jordan Mirocha <mirochaj@gmail.com>
"""
try:
    import atpy
    atpyOK = True
except ImportError:
    atpyOK = False
from .. import units
import numpy as np

def open_1d_txt(filename, xaxcol=0, datacol=1, errorcol=2, **kwargs):
    """
    Attempt to read a 1D spectrum from a text file assuming
    wavelength as the first column, data as the second, and
    (optionally) error as the third.  

    kwargs are passed to atpy.Table
    """
    if not atpyOK:
        print "WARNING: atpy not installed; will use simple reader instead."
        
        data, error, XAxis, T = simple_txt(filename, xaxcol = xaxcol, 
            datacol = datacol, errorcol = errorcol, **kwargs)
        
    else:
        T = atpy.Table(filename, type='ascii', masked=True, **kwargs)
        
        xarr = T.data[T.data.dtype.names[xaxcol]]
        data = T.data[T.data.dtype.names[datacol]]
        if len(T.columns) > errorcol:
            error = T.data[T.data.dtype.names[errorcol]]
        else:
            # assume uniform, zero error
            error = data*0 

        if 'xunits' in T.keywords:
            xunits = T.keywords['xunits']
        else:
            xunits = 'unknown'

        XAxis = units.SpectroscopicAxis(xarr,xunits)
    
    # Need this in Spectrum class to correctly parse header    
    T.xaxcol = xaxcol
    T.datacol = datacol 

    return data, error, XAxis, T

def simple_txt(filename, xaxcol=0, datacol=1, errorcol=2, **kwargs):
    """
    Very simple method for reading columns from ASCII file.
    """
    
    f = open(filename, 'r')
    
    colunits = []
    coldata = []
    for i, line in enumerate(f):
        
        # Ignore blank lines
        if not line.strip(): continue
        
        # Possibly read in header
        if line.split()[0][0] == '#': 
            if i == 0: hdr = line.split()[1:]
            if i == 1: colunits = line.split()[1:]
            
            continue
        
        coldata.append(line.split())
        
        for j, element in enumerate(coldata[-1]):
            try: coldata[-1][j] = float(element)
            except ValueError: coldata[-1][j] = str(element)
                
    f.close()
    
    if not colunits:
        colunits = ['unknown', 'unknown', 'unknown']
    if not hdr:
        hdr = ['unknown', 'unknown', 'unknown']    
        
    N = len(hdr)
    coldata = zip(*coldata)     
            
    # Prepare to return data
    data = coldata[datacol]
    xarr = coldata[xaxcol]
    if errorcol > N - 1:
        error = data*0
    else:
        error = coldata[errorcol]

    XAxis = units.SpectroscopicAxis(xarr, colunits[xaxcol]) 
    
    # Create atPy style Table instance
    T = dummy_class()
    T.data = dummy_class()
    T.data.dtype = dummy_class()
    T.data.dtype.names = hdr
    T.columns = {}
    T.columns[T.data.dtype.names[xaxcol]] = dummy_class()
    T.columns[T.data.dtype.names[xaxcol]].unit = colunits[xaxcol]
    T.columns[T.data.dtype.names[datacol]] = dummy_class()
    T.columns[T.data.dtype.names[datacol]].unit = colunits[datacol]
                
    return np.array(data), np.array(error), XAxis, T
    
class dummy_class:
    def __init__(self):
        pass


