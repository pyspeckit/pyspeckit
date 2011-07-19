import pyfits
from .. import units
import numpy as np
import numpy.ma as ma

def tspec_reader(filename):
    """ 
    Read TripleSpec reduced data

    It may be reduced by either SPEXTOOL or IRAF.

    AFAIK, Spextool will use the INSTR keyword, IRAF will use INSTRUME
    """
    fitsfile = pyfits.open(filename)
    header = fitsfile[0].header
    data = fitsfile[0].data

    if header.get('INSTR') == 'APO Triplespec':
        # read in SPEXTOOL spectrum
        xarr = np.array(data[0,:],dtype='float64')
        spec = ma.array(data[1,:],dtype='float64')
        spec.mask = spec!=spec
        errspec = ma.array(data[2,:],dtype='float64')
        errspec.mask = (spec!=spec)+(errspec!=errspec)
        xunits = header.get('XUNITS')
        unit = header.get('YUNITS')
        header.update('BUNIT',unit)
        xtype = 'wavelength'
    elif header.get('INSTRUME') == 'tspec':
        dv,v0,p3 = header['CD1_1'],header['CRVAL1'],header['CRPIX1']
        header.update('CDELT1',dv)
        xconv = lambda v: ((v-p3+1)*dv+v0)
        xarr = xconv(np.arange(len(spec)))
        wat = dict([s.split("=") for s in header.get('WAT1_001').split()])
        xunits = wat['units']
        xtype = wat['label']
        #unit = header.get('BUNIT')
        
    XAxis = units.SpectroscopicAxis(xarr,xunits,xtype=xtype)

    return spec,errspec,XAxis,header
