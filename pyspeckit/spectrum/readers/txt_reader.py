try:
    import atpy
    atpyOK = True
except ImportError:
    atpyOK = False
from .. import units

def open_1d_txt(filename,errspecnum=2):
    if not atpyOK: print "atpy not installed; cannot read txt files"
    else:
        T = atpy.Table(filename, type='ascii',
                Reader=atpy.asciitables.asciitable.CommentedHeader, masked=True)

        
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

