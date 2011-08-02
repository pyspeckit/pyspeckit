try:
    import atpy
    atpyOK = True
except ImportError:
    atpyOK = False
from .. import units

def open_1d_txt(filename,errspecnum=2, **kwargs):
    """
    Attempt to read a 1D spectrum from a text file assuming
    wavelength as the first column, data as the second, and
    (optionally) error as the third.  

    kwargs are passed to atpy.Table
    """
    if not atpyOK:
        print "atpy not installed; cannot read txt files"
    else:
        T = atpy.Table(filename, type='ascii', masked=True, **kwargs)
        
        xarr = T.data[T.data.dtype.names[0]]
        data = T.data[T.data.dtype.names[1]]
        if len(T.columns) > errspecnum:
            error = T.data[T.data.dtype.names[errspecnum]]
        else:
            # assume uniform, nonzero error
            error = data*0 + 1.0

        if 'xunits' in T.keywords:
            xunits = T.keywords['xunits']
        else:
            xunits = 'unknown'

        XAxis = units.SpectroscopicAxis(xarr,xunits)

        return data,error,XAxis,T

