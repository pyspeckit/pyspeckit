from __future__ import print_function
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
from .. import units
import numpy as np
import numpy.ma as ma

def tspec_reader(filename, merged=True, specnum=0, **kwargs):
    """ 
    Read TripleSpec or SPEXTOOL reduced data

    It may be reduced by either SPEXTOOL or IRAF.

    AFAIK, Spextool will use the INSTR keyword, IRAF will use INSTRUME

    non-merged SPEXTOOL files are acceptable; they'll probably be read correctly but no guarantees

    kwargs is ignored
    """
    fitsfile = pyfits.open(filename)
    header = fitsfile[0].header
    data = fitsfile[0].data

    if merged:
        if header.get('INSTR') == 'APO Triplespec' or header.get('INSTR') == 'SpeX':
            # read in SPEXTOOL spectrum
            xarr = np.array(data[0,:],dtype='float64')
            spec = ma.array(data[1,:],dtype='float64')
            spec.mask = spec!=spec
            errspec = ma.array(data[2,:],dtype='float64')
            errspec.mask = (spec!=spec)+(errspec!=errspec)
            xunits = header.get('XUNITS')
            unit = header.get('YUNITS')
            header['BUNIT'] = unit
            xtype = 'wavelength'
        elif header.get('INSTRUME') == 'tspec':
            dv,v0,p3 = header['CD1_1'],header['CRVAL1'],header['CRPIX1']
            header['CDELT1'] = dv
            xconv = lambda v: ((v-p3+1)*dv+v0)
            spec = ma.array(data[1,:],dtype='float64')
            xarr = xconv(np.arange(len(spec)))
            wat = dict([s.split("=") for s in header.get('WAT1_001').split()])
            xunits = wat['units']
            xtype = wat['label']
            #unit = header.get('BUNIT')
        else: # try this...
            xarr = np.array(data[specnum,0,:],dtype='float64')
            spec = ma.array(data[specnum,1,:],dtype='float64')
            spec.mask = spec!=spec
            errspec = ma.array(data[specnum,2,:],dtype='float64')
            errspec.mask = (spec!=spec)+(errspec!=errspec)
            xunits = header.get('XUNITS')
            unit = header.get('YUNITS')
            header['BUNIT'] = unit
            xtype = 'wavelength'
    else:
        xarr = np.array(data[specnum,0,:],dtype='float64')
        spec = ma.array(data[specnum,1,:],dtype='float64')
        spec.mask = spec!=spec
        errspec = ma.array(data[specnum,2,:],dtype='float64')
        errspec.mask = (spec!=spec)+(errspec!=errspec)
        xunits = header.get('XUNITS')
        unit = header.get('YUNITS')
        header['BUNIT'] = unit
        xtype = 'wavelength'
        
    XAxis = units.SpectroscopicAxis(xarr,xunits,xtype=xtype)

    return spec,errspec,XAxis,header
