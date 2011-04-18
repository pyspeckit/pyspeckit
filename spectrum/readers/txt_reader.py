try:
    import atpy
    atpyOK = True
except ImportError:
    atpyOK = False
import spectrum.units as units

def open_1d_txt(filename):
    if not atpyOK: print "atpy not installed; cannot read txt files"
    else:
        T = atpy.Table(filename, type='ascii',
                Reader=atpy.asciitables.asciitable.CommentedHeader, masked=True)

        
        xarr = T.data[T.data.dtype.names[0]]
        data = T.data[T.data.dtype.names[1]]
        if len(T.columns) > 2:
            error = T.data[T.data.dtype.names[2]]
        else:
            # assume uniform, nonzero error
            error = data*0 + 1.0

        if 'xunits' in T.keywords:
            xunits = T.keywords['xunits']
        else:
            xunits = 'unknown'

        XAxis = units.SpectroscopicAxis(xarr,xunits)

        return data,error,XAxis,T

