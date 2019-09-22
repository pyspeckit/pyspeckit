from __future__ import print_function
from six import iteritems
import numpy as np
import warnings
r2deg = 180/np.pi
"""
.. TODO::
    When any section length is zero, that means the following values are to be
    ignored.  No warning is needed.
"""

def read_lmv(fn):
    """
    Read an LMV cube file

    Specification is primarily in GILDAS image_def.f90
    """

    header = {}
    with open(fn,'rb') as lf:
        # lf for "LMV File"
        filetype = lf.read(12)
        if filetype != 'GILDAS-IMAGE':
            raise TypeError("File is not a GILDAS Image file")

        # fmt probably matters!  Default is "r4", i.e. float32 data, but could be float64
        fmt = np.fromfile(lf, dtype='int32', count=1) # 4
        # number of data blocks
        ndb = np.fromfile(lf, dtype='int32', count=1) # 5
        gdf_type = np.fromfile(lf, dtype='int32', count=1) # 6
        # Reserved Space
        reserved_fill = np.fromfile(lf, dtype='int32', count=4) # 7
        general_section_length = np.fromfile(lf, dtype='int32', count=1) # 11
        #print "Format: ",fmt," ndb: ",ndb, " fill: ",fill," other: ",unknown

        # pos 12
        naxis,naxis1,naxis2,naxis3,naxis4 = np.fromfile(lf,count=5,dtype='int32')
        header['NAXIS'] = naxis
        header['NAXIS1'] = naxis1
        header['NAXIS2'] = naxis2
        header['NAXIS3'] = naxis3
        header['NAXIS4'] = naxis4

        # We are indexing bytes from here; CLASS indices are higher by 12
        # pos 17
        header['CRPIX1'] = np.fromfile(lf,count=1,dtype='float64')[0]
        header['CROTA1'] = np.fromfile(lf,count=1,dtype='float64')[0]
        header['CRVAL1'] = np.fromfile(lf,count=1,dtype='float64')[0] * r2deg
        header['CRPIX2'] = np.fromfile(lf,count=1,dtype='float64')[0]
        header['CROTA2'] = np.fromfile(lf,count=1,dtype='float64')[0]
        header['CRVAL2'] = np.fromfile(lf,count=1,dtype='float64')[0] * r2deg
        header['CRPIX3'] = np.fromfile(lf,count=1,dtype='float64')[0]
        header['CROTA3'] = np.fromfile(lf,count=1,dtype='float64')[0]
        header['CRVAL3'] = np.fromfile(lf,count=1,dtype='float64')[0]
        header['CRPIX4'] = np.fromfile(lf,count=1,dtype='float64')[0]
        header['CROTA4'] = np.fromfile(lf,count=1,dtype='float64')[0]
        header['CRVAL4'] = np.fromfile(lf,count=1,dtype='float64')[0]
        # pos 41
        #print "Post-crval",lf.tell()
        blank_section_length = np.fromfile(lf,count=1,dtype='int32')
        if blank_section_length != 8:
            warnings.warn("Invalid section length found for blanking section")
        header['BLANK'] = np.fromfile(lf,count=1,dtype='float32')[0] # 42
        header['TOLERANC'] = np.fromfile(lf,count=1,dtype='int32')[0] # 43 eval = tolerance
        extrema_section_length = np.fromfile(lf,count=1,dtype='int32')[0] # 44
        if extrema_section_length != 40:
            warnings.warn("Invalid section length found for extrema section")
        vmin,vmax = np.fromfile(lf,count=2,dtype='float32') # 45
        xmin,xmax,ymin,ymax,zmin,zmax = np.fromfile(lf,count=6,dtype='int32') # 47
        wmin,wmax = np.fromfile(lf,count=2,dtype='int32') # 53
        description_section_length = np.fromfile(lf,count=1,dtype='int32')[0] # 55
        if description_section_length != 72:
            warnings.warn("Invalid section length found for description section")
        #strings = lf.read(description_section_length) # 56
        header['BUNIT'] = lf.read(12) # 56
        header['CTYPE1'] = lf.read(12) # 59
        header['CTYPE2'] = lf.read(12) # 62
        header['CTYPE3'] = lf.read(12) # 65
        header['CTYPE4'] = lf.read(12) # 68
        header['COOSYS'] = lf.read(12) # 71
        position_section_length = np.fromfile(lf,count=1,dtype='int32') # 74
        if position_section_length != 48:
            warnings.warn("Invalid section length found for position section")
        header['OBJNAME'] = lf.read(4*3) # 75
        header['RA'] = np.fromfile(lf, count=1, dtype='float64')[0] * r2deg # 78
        header['DEC'] = np.fromfile(lf, count=1, dtype='float64')[0] * r2deg # 80
        header['GLON'] = np.fromfile(lf, count=1, dtype='float64')[0] * r2deg # 82
        header['GLAT'] = np.fromfile(lf, count=1, dtype='float64')[0] * r2deg # 84
        header['EQUINOX'] = np.fromfile(lf,count=1,dtype='float32')[0] # 86
        header['PROJWORD'] = lf.read(4) # 87
        header['PTYP'] = np.fromfile(lf,count=1,dtype='int32')[0] # 88
        header['A0'] = np.fromfile(lf,count=1,dtype='float64')[0] # 89
        header['D0'] = np.fromfile(lf,count=1,dtype='float64')[0] # 91
        header['PANG'] = np.fromfile(lf,count=1,dtype='float64')[0] # 93
        header['XAXI'] = np.fromfile(lf,count=1,dtype='float32')[0] # 95
        header['YAXI'] = np.fromfile(lf,count=1,dtype='float32')[0] # 96
        spectroscopy_section_length = np.fromfile(lf,count=1,dtype='int32') # 97
        if spectroscopy_section_length != 48:
            warnings.warn("Invalid section length found for spectroscopy section")
        header['RECVR'] = lf.read(12) # 98
        header['FRES'] = np.fromfile(lf,count=1,dtype='float64')[0] # 101
        header['IMAGFREQ'] = np.fromfile(lf,count=1,dtype='float64')[0] # 103 "FIMA"
        header['REFFREQ'] = np.fromfile(lf,count=1,dtype='float64')[0] # 105
        header['VRES'] = np.fromfile(lf,count=1,dtype='float32')[0] # 107
        header['VOFF'] = np.fromfile(lf,count=1,dtype='float32')[0] # 108
        header['FAXI'] = np.fromfile(lf,count=1,dtype='int32')[0] # 109
        resolution_section_length = np.fromfile(lf,count=1,dtype='int32')[0] # 110
        if resolution_section_length != 12:
            warnings.warn("Invalid section length found for resolution section")
        #header['DOPP'] = np.fromfile(lf,count=1,dtype='float16')[0] # 110a ???
        #header['VTYP'] = np.fromfile(lf,count=1,dtype='int16')[0] # 110b
        # integer, parameter :: vel_unk = 0      ! Unsupported referential :: planetary...)
        # integer, parameter :: vel_lsr = 1      ! LSR referential
        # integer, parameter :: vel_hel = 2      ! Heliocentric referential
        # integer, parameter :: vel_obs = 3      ! Observatory referential
        # integer, parameter :: vel_ear = 4      ! Earth-Moon barycenter referential
        # integer, parameter :: vel_aut = -1     ! Take referential from data
        header['BMAJ'] = np.fromfile(lf,count=1,dtype='float32')[0] # 111
        header['BMIN'] = np.fromfile(lf,count=1,dtype='float32')[0] # 112
        header['BPA'] = np.fromfile(lf,count=1,dtype='float32')[0] # 113
        noise_section_length = np.fromfile(lf,count=1,dtype='int32')
        if noise_section_length != 0:
            warnings.warn("Invalid section length found for noise section")
        header['NOISE'] = np.fromfile(lf,count=1,dtype='float32')[0] # 115
        header['RMS'] = np.fromfile(lf,count=1,dtype='float32')[0] # 116
        astrometry_section_length = np.fromfile(lf,count=1,dtype='int32')
        if astrometry_section_length != 0:
            warnings.warn("Invalid section length found for astrometry section")
        header['MURA'] = np.fromfile(lf,count=1,dtype='float32')[0] # 118
        header['MUDEC'] = np.fromfile(lf,count=1,dtype='float32')[0] # 119
        header['PARALLAX'] = np.fromfile(lf,count=1,dtype='float32')[0] # 120

        other_info = np.fromfile(lf, count=7, dtype='float32') # 121-end
        if not np.all(other_info == 0):
            warnings.warn("Found additional information in the last 7 bytes")

        endpoint = 508
        if lf.tell() != endpoint:
            raise ValueError("Header was not parsed correctly")

        data = np.fromfile(lf, count=naxis1*naxis2*naxis3, dtype='float32')

    data[data == header['BLANK']] = np.nan

    # for no apparent reason, y and z are 1-indexed and x is zero-indexed
    if (wmin-1,zmin-1,ymin-1,xmin) != np.unravel_index(np.nanargmin(data),
                                                       [naxis4,naxis3,naxis2,naxis1]):
        warnings.warn("Data min location does not match that on file.  "
                      "Possible error reading data.")
    if (wmax-1,zmax-1,ymax-1,xmax) != np.unravel_index(np.nanargmax(data),
                                                       [naxis4,naxis3,naxis2,naxis1]):
        warnings.warn("Data max location does not match that on file.  "
                      "Possible error reading data.")
    if np.nanmax(data) != vmax:
        warnings.warn("Data max does not match that on file.  "
                      "Possible error reading data.")
    if np.nanmin(data) != vmin:
        warnings.warn("Data min does not match that on file.  "
                      "Possible error reading data.")

    return data.reshape([naxis4,naxis3,naxis2,naxis1]),header
    # debug
    #return data.reshape([naxis3,naxis2,naxis1]), header, hdr_f, hdr_s, hdr_i, hdr_d, hdr_d_2

def read_lmv_tofits(fn):
    from astropy.io import fits
    data,header = read_lmv(fn)

    cards = [fits.header.Card(k,v) for k,v in iteritems(header)]
    Header = fits.Header(cards)
    hdu = fits.PrimaryHDU(data=data, header=Header)
    return hdu
