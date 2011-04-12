import pyfits
import numpy
import struct
from numpy import pi
"""
A note from the CLASS developers: 

> Dear Adam,
> 
> Our policy regarding the Class binary data format changed 5 years again. It now
> is private to class developers. We partly document it for the sake of
> transparency (and it happens that the current documentation is out of date
> because of lack of manpower) but only the API is public. The reason is that it
> leaves us the freedom (even though seldom used) to make it evolve with the new
> needs of the instruments. In other words, there is no warranty that the class
> binary data format will stay as it is today but we always ensure that the Class
> library is able to read a Class file, whatever its age: backward compatibility
> is a must-have feature for us.
> 
> Our recommendation is thus that you use by a mean or an other the Class
> library. Since you are a Python user, it should be easy for you to use the
> Gildas-Python binding. You will find at the end of this email an example of a
> session which uses it. We hope this will help you.
> 
> Best regards,
> 
> Sebastien Bardeau and Jerome Pety

Because I want purely non-proprietary data access, I decided to go ahead and
try to figure out what a CLASS file consists of anyway.

The structure is pretty complicated - each observation has an INDEX and a
HEADER.  The INDEX is what CLASS's LIST command accesses - it has scan numbers,
source numbers, and "telescope" names (the quotes are because "telescope" often
refers to the instrument on the telescope, not the dish).

The various header sections include very useful information, but header section
-2 doesn't appear to exist in the data (Even though it is stated to be in the data
at the tops of the headers).  The sort of stuff I haven't been able to parse looks like
SHM_MAP: Segment key 0xcafe015, mapped to shmid 0x44e48010o shm 0xcafe015g 2, for spec 0 mode        Sky bin 0 walsh

Without testing, I'm not sure about certain parameters, such as lamof/betof -
are they in radians (like lam/bet) or some other system?

Additionally, the projection system is unclear to me... there is a 'systex'
array in titout.f90 that has 6 entries, but the index in the smt files is 7.
So...  I guess I have to assume simple fk5, tangent?

"""

""" Specification: http://iram.fr/IRAMFR/GILDAS/doc/html/class-html/node58.html """
filetype_dict = {'1A  ':'Multiple_IEEE','1   ':'Multiple_Vax','1B  ':'Multiple_EEEI',
    '9A  ':'Single_IEEE','9   ':'Single_Vax','9B  ':'Single_EEEI'}

keys_lengths = { -2: [
     ( 'NUM'     ,1,'int32'), # Observation number
     ( 'VER'     ,1,'int32'), # Version number
     ( 'TELES'   ,3,'|S12') , # Telescope name
     ( 'DOBS'    ,1,'int32'), # Date of observation
     ( 'DRED'    ,1,'int32'), # Date of reduction
     ( 'TYPEC'   ,1,'int32'), # Type of coordinates
     ( 'KIND'    ,1,'int32'), # Type of data
     ( 'QUAL'    ,1,'int32'), # Quality of data
     ( 'SCAN'    ,1,'int32'), # Scan number
     ( 'SUBSCAN' ,1,'int32'), # Subscan number
     ( 'UT'      ,2,'float64'), #  rad UT of observation
     ( 'ST'      ,2,'float64'), #  rad LST of observation
     ( 'AZ'      ,1,'float32'), #  rad Azimuth
     ( 'EL'      ,1,'float32'), #  rad Elevation
     ( 'TAU'     ,1,'float32'), # neper Opacity
     ( 'TSYS'    ,1,'float32'), #    K System temperature
     ( 'TIME'    ,1,'float32'), #    s Integration time
     ( 'XUNIT'   ,1,'int32'),   # code X unit (if xcoord_sec is present)
     ] ,
     'POSITION': [
    ('SOURC',3,'|S12')  , #  [ ] Source name
    ('EPOCH',1,'float32'), #  [ ] Epoch of coordinates
    ('LAM'  ,2,'float64'), #[rad] Lambda
    ('BET'  ,2,'float64'), #[rad] Beta
    ('LAMOF',1,'float32'), #  [rad] Offset in Lambda
    ('BETOF',1,'float32'), #  [rad] Offset in Beta
    ('PROJ' ,1,'int32')  , # [rad] Projection system
    ],
     'SPECTRO': [
     #('align'  ,1,'int32'),   #  [    ] Alignment padding
     ('LINE'   ,3,'|S12'),    #  [    ] Line name
     ('RESTF'  ,2,'float64'), #  [ MHz] Rest frequency
     ('NCHAN'  ,1,'int32'),   #  [    ] Number of channels
     ('RCHAN'  ,1,'float32'), #  [    ] Reference channels
     ('FRES'   ,1,'float32'), #  [ MHz] Frequency resolution
     ('FOFF'   ,1,'float32'), #  [ MHz] Frequency offset
     ('VRES'   ,1,'float32'), #  [km/s] Velocity resolution
     ('VOFF'   ,1,'float32'), #  [km/s] Velocity at reference channel
     ('BAD'    ,1,'float32'), #  [    ] Blanking value
     #('ALIGN_1',1,'int32'),   #  [    ] Alignment padding
     ('IMAGE'  ,2,'float64'), #  [ MHz] Image frequency
     #('ALIGN_2',1,'int32'),   #  [    ] Alignment padding
     ('VTYPE'  ,1,'int32'),   #  [code] Type of velocity
     ('DOPPLER',2,'float64'), #  [    ] Doppler factor = -V/c (CLASS convention)
     ],
     'CALIBRATION': [
     ('ALIGN',1,'int32'),    # BUFFER
     ('BEEFF',1,'float32'),   # [ ] Beam efficiency
     ('FOEFF',1,'float32'),   # [ ] Forward efficiency
     ('GAINI',1,'float32'),   # [ ] Image/Signal gain ratio
     ('H2OMM',1,'float32'),   # [ mm] Water vapor content
     ('PAMB',1,'float32'),   # [ hPa] Ambient pressure
     ('TAMB',1,'float32'),   # [ K] Ambient temperature
     ('TATMS',1,'float32'),   # [ K] Atmosphere temp. in signal band
     ('TCHOP',1,'float32'),   # [ K] Chopper temperature
     ('TCOLD',1,'float32'),   # [ K] Cold load temperature
     ('TAUS',1,'float32'),   # [neper] Opacity in signal band
     ('TAUI',1,'float32'),   # [neper] Opacity in image band
     ('TATMI',1,'float32'),   # [ K] Atmosphere temp. in image band
     ('TREC',1,'float32'),   # [ K] Receiver temperature
     ('CMODE',1,'int32'),   # [ code] Calibration mode
     ('ATFAC',1,'float32'),   # [ ] Applied calibration factor
     ('ALTI',1,'float32'),   # [ m] Site elevation
     ('COUNT(3)',3,'3float32'),   # [count] Power of Atm., Chopp., Cold
     ('LCALOF',1,'float32'),   # [ rad] Longitude offset for sky measurement
     ('BCALOF',1,'float32'),   # [ rad] Latitude offset for sky measurement
     #('GEOLONG',2,'float64'),   # [ rad] Geographic longitude of observatory
     #('GEOLAT',2,'float64'),   # [ rad] Geographic latitude of observatory
         ],
     
     }

def _read_bytes(f, n):
    '''Read the next `n` bytes (from idlsave)'''
    return f.read(n)


def _read_byte(f):
    '''Read a single byte (from idlsave)'''
    return numpy.uint8(struct.unpack('>B', f.read(4)[:1])[0])

def _read_int16(f):
    '''Read a signed 16-bit integer (from idlsave)'''
    return numpy.int16(struct.unpack('>h', f.read(4)[2:4])[0])

def _read_int32(f):
    '''Read a signed 32-bit integer (from idlsave)'''
    return numpy.int32(struct.unpack('>i', f.read(4))[0])

def _read_float32(f):
    '''Read a 32-bit float (from idlsave)'''
    return numpy.float32(struct.unpack('>f', f.read(4))[0])

def _align_32(f):
    '''Align to the next 32-bit position in a file (from idlsave)'''

    pos = f.tell()
    if pos % 4 != 0:
        f.seek(pos + 4 - pos % 4)
    return

def read_word(f,length):
    if length > 0:
        chars = _read_bytes(f, length)
        _align_32(f)
    else:
        chars = None
    return chars

def read_int(f):
    return struct.unpack('i',f.read(4))

def read_index(f):
    index = {
                "XBLOC":_read_byte(f),
                "XNUM":_read_byte(f),
                "XVER":_read_byte(f),
                "XSOURC":read_word(f,12),
                "XLINE":read_word(f,12),
                "XTEL":read_word(f,12),
                "XDOBS":_read_byte(f),
                "XDRED":_read_byte(f),
                "XOFF1":_read_float32(f),# 		 first offset (real, radians) 
                "XOFF2":_read_float32(f),# 		 second offset (real, radians) 
                "XTYPE":read_word(f,2),# 		 coordinate system ('EQ'', 'GA', 'HO') 
                "XKIND":_read_byte(f),# 		 Kind of observation (0: spectral, 1: continuum, ) 
                "XQUAL":_read_byte(f),# 		 Quality (0-9)  
                "XSCAN":_read_byte(f),# 		 Scan number 
                "XPOSA":_read_byte(f),# 		 Position Angle 
                "XFRONT":read_word(f,8),# 		 (8 char) Front-end  ID (PROPOSED) 
                "XBACK" :read_word(f,8),# 		 (8 char) Back-end   ID (PROPOSED) 
                "XPROC" :read_word(f,8),# 		 (8 char) Procedure  ID (PROPOSED) 
                "XPROJ" :read_word(f,8),# 		 (8 char) Project    ID (PROPOSED) 
                "UNUSED":read_word(f,6),
                "BLANK" :read_word(f,4)
            }
    return index

def read_header(f,type=0):
    if type in keys_lengths:
        hdrsec = [(x[0],numpy.fromfile(f,count=1,dtype=x[2])[0])
                for x in keys_lengths[type]]
        return dict(hdrsec)
    pass


def read_obshead(f,verbose=False):
    #IDcode = f.read(4)
    #if verbose: print "IDcode: ",IDcode
    #if IDcode.strip() != '2':
    #    raise IndexError("Failure at %i" % (f.tell() - 4))
    nblocks,nbyteob,data_address,nheaders,data_length,obindex,nsec,obsnum= numpy.fromfile(f,count=8,dtype='int32')
    if verbose: print "nblocks,nbyteob,data_address,data_length,nheaders,obindex,nsec,obsnum",nblocks,nbyteob,data_address,data_length,nheaders,obindex,nsec,obsnum
    seccodes = numpy.fromfile(f,count=nsec,dtype='int32')
    secaddr = numpy.fromfile(f,count=nsec,dtype='int32')
    seclen = numpy.fromfile(f,count=nsec,dtype='int32')
    if verbose: print "Section codes, addresses, lengths: ",seccodes,secaddr,seclen

    hdr = {'NBLOCKS':nblocks, 'NBYTEOB':nbyteob, 'DATA_ADDRESS':data_address,
            'DATA_LENGTH':data_length, 'NHEADERS':nheaders, 'OBINDEX':obindex,
            'NSEC':nsec, 'OBSNUM':obsnum, 'SECCODES':seccodes,
            'SECADDR':secaddr, 'SECLEN':seclen}

    #return obsnum,seccodes
    return obsnum,hdr


def read_class(filename, DEBUG=False):
    f = open(filename,'rb')
    filelen = len(f.read())
    f.seek(0)

    filetype = f.read(4)
    if filetype in filetype_dict:
        print "File is type %s" % filetype_dict[filetype]
    else:
        raise TypeError("File type error: %s." % filetype)
    
    nextblock,nindex,nex,nrecords = numpy.fromfile(f,count=4,dtype='int32')
    if DEBUG: print "nextblock,nindex,nex,nrecords",nextblock,nindex,nex,nrecords
    firstblock=numpy.fromfile(f,count=nex,dtype='int32') #[_read_byte(f) for ii in xrange(nex)]

    if DEBUG: print "firstblock",firstblock
    if DEBUG: print "Done with header stuff at position %i" % f.tell()

    indexes = []
    for ii in xrange(nindex):
        f.seek(128*(ii)+(firstblock[0]-1)*512)
        index = read_index(f)
        # OLD DEBUG if index['XLINE'] not in ('HCOP(3-2)   ','N2HP(3-2)   '):
        # OLD DEBUG     raise Exception("Stopped at %i" % ii)
        indexes.append(index)

    #if f.tell() % 128 != 0:
    #    f.seek((f.tell()/128 + 1)*128)

    spectra = []
    header  = []
    spcount = 0
    jj = -1
    while f.tell() < filelen:
        jj += 1
        startpos = f.tell()
        IDcode = f.read(4)
        if IDcode.strip() != '2':
            f.seek(startpos)
            index = read_index(f)
            if 'HCO' in index['XLINE'] or 'N2H' in index['XLINE']:
                if DEBUG: print "Found an index at %i: " % startpos,index
                indexes.append(index)
                continue
            else:
                print "Failure at %i: %i (/512=%i)" % (jj,f.tell(),f.tell()/512)
                raise Exception("Failure at %i: %i" % (jj,f.tell()))
        else:
            if DEBUG: print "Reading HEADER and OBSERVATION at %i" % f.tell()
            obsnum,obshead = read_obshead(f,verbose=DEBUG)
        pos = f.tell()
        if DEBUG:
            print " UNCLEAR INFO AT %i PAST %i:" % (f.tell()-startpos,startpos),f.read(168) 
            f.seek(pos)
            #f.seek(startpos)
            #print numpy.fromfile(f,count=168/4+10,dtype='int32') 
            #f.seek(startpos)
            #print numpy.fromfile(f,count=168/4+10,dtype='float32') 
            #f.seek(startpos)
            #print numpy.fromfile(f,count=168+40,dtype='int8') 
        Header2 = read_header(f,type=-2)
        if f.tell() != pos + 168:
            #print "Wrong position %i, skipping to %i" % (f.tell(),pos+168)
            f.seek(pos+168)
        Header3 = read_header(f,type='POSITION')
        Header4 = read_header(f,type='SPECTRO')
        if DEBUG: print "Line %i (byte %i) - OBSERVATION %i (%i): %s, %s" % (f.tell(),(f.tell()-startpos)/4,obsnum,spcount,Header3['SOURC'],Header4['LINE']),
        Header14 = read_header(f,type='CALIBRATION')
        hdr = Header2
        hdr.update(Header3)
        hdr.update(Header4)
        hdr.update(Header14)
        #hdr.update(unclear)
        hdr.update(obshead)
        hdr.update({'OBSNUM':obsnum,'RECNUM':spcount})
        hdr.update({'RA':hdr['LAM']/pi*180,'DEC':hdr['BET']/pi*180})
        hdr.update({'RAoff':hdr['LAMOF']/pi*180,'DECoff':hdr['BETOF']/pi*180})
        hdr.update({'OBJECT':hdr['SOURC'].strip()})
        hdr.update({'BUNIT':'Tastar'})
        hdr.update({'EXPOSURE':hdr['TIME']})
        spcount += 1

        nchan = hdr['NCHAN']
        if nchan != hdr['DATA_LENGTH']:
            raise ValueError("data_length != nchan")
        if DEBUG: print "Spectrum has %i channels at %i" % (nchan,f.tell())
        spectrum = numpy.fromfile(f,count=nchan,dtype='float32')
        if DEBUG > 2: print "First digits of spectrum: ",spectrum[:10]
        spectra.append( spectrum )
        header.append( hdr )
        f.seek((f.tell()/nchan + 1)*nchan)

    f.close()
    return spectra,header,indexes

import spectrum
def make_axis(header):

    rest_frequency = header.get('RESTF')
    xunits = 'MHz'
    nchan = header.get('NCHAN')
    voff = header.get('VOFF')
    foff = header.get('FOFF')
    doppler = header.get('DOPPLER')
    fres = header.get('FRES')
    refchan = header.get('RCHAN')
    imfreq = header.get('IMAGE')

    xarr = (numpy.arange(nchan) - refchan) * fres + rest_frequency
    XAxis = spectrum.units.SpectroscopicAxis(xarr,'MHz',frame='rest',reffreq=rest_frequency)

    return XAxis
    
def class_to_obsblocks(filename,telescope,line):
    """
    Load an entire CLASS observing session into a list of ObsBlocks based on
    matches to the 'telescope' and 'line' names
    """
    spectra,header,indexes = read_class(filename)

    H = pyfits.Header()

    obslist = []
    lastscannum = -1
    spectrumlist = None
    for sp,hdr,ind in zip(spectra,header,indexes):
        hdr.update(ind)
        # this is slow but necessary...
        for k,v in hdr.iteritems():
            try:
                H.update(k,v)
            except ValueError:
                pass
                #H.update(k,str(v))
        scannum = hdr['XSCAN']
        if hdr['XTEL'].strip() not in telescope:
            continue
        if hdr['LINE'].strip() not in line:
            continue

        if scannum != lastscannum:
            lastscannum = scannum
            if spectrumlist is not None:
                obslist.append(spectrum.ObsBlock(spectrumlist))
            xarr = make_axis(hdr)
            spectrumlist = [(
                spectrum.Spectrum(xarr=xarr,
                    header=H,
                    data=sp))]
        else:
            spectrumlist.append(
                spectrum.Spectrum(xarr=xarr,
                    header=H,
                    data=sp))

    return obslist

def class_to_spectra(filename):
    """
    Load each individual spectrum within a CLASS file into a list of Spectrum
    objects
    """
    spectra,header,indexes = read_class(filename)

    spectrumlist = []
    for sp,hdr,ind in zip(spectra,header,indexes):
        hdr.update(ind)
        xarr = make_axis(hdr)
        spectrumlist.append(
            spectrum.Spectrum(xarr=xarr,
                header=hdr,
                data=sp))

    return spectrumlist

if __name__ == "__main__":
    fn1 = '/Users/adam/work/bolocam/hht/class_001.smt' 
    #fn1 = '/Users/adam/work/bolocam/hht/test_SMT-F1M-VU-20824-073.cls' 
    #fn2 = '/Users/adam/work/bolocam/hht/test_SMT-F1M-VU-79472+203.cls'
    F1 = read_class(fn1)#,DEBUG=True)
    #F2 = read_class(fn2)
    n2hp = class_to_obsblocks(fn1,telescope=['SMT-F1M-HU','SMT-F1M-VU'],line='N2HP(3-2)')
    hcop = class_to_obsblocks(fn1,telescope=['SMT-F1M-HL','SMT-F1M-VL'],line='HCOP(3-2)')
