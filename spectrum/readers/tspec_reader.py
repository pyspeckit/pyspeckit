import pyfits
import spectrum.units as units
import numpy as np

def tspec_reader(filename):
    """ 
    Read TripleSpec reduced data

    It may be reduced by either SPEXTOOL or IRAF.

    AFAIK, Spextool will use the INSTR keyword, IRAF will use INSTRUME
    """
    fitsfile = pysits.open(filename)
    header = fitsfile[0].header
    data = fitsfile[0].data

    if header.get('INSTR') == 'APO Triplespec':
        # read in SPEXTOOL spectrum
        xarr = data[0,:]
        spec = data[1,:]
        errspec = data[2,:]
        xunits = header.get('XUNITS')
        units = header.get('YUNITS')
        header.update('BUNIT',units)
        xtype = 'wavelength'
    elif header.get('INSTRUME') == 'tspec':
        dv,v0,p3 = hdr['CD1_1'],hdr['CRVAL1'],hdr['CRPIX1']
        hdr.update('CDELT1',dv)
        xconv = lambda v: ((v-p3+1)*dv+v0)
        xarr = xconv(np.arange(len(spec)))
        wat = dict([s.split("=") for s in header.get('WAT1_001').split()])
        xunits = wat['units']
        xtype = wat['label']
        units = header.get('BUNIT')
        
    XAxis = units.SpectroscopicAxis(xarr,xunits,xtype=xtype,reffreq=reffreq)

    return spec,errspec,XAxis,hdr
