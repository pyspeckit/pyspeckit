"""
------------------------
GILDAS CLASS file reader
------------------------

Read a CLASS file into an :class:`pyspeckit.spectrum.ObsBlock`
"""
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
import numpy
import numpy as np
from numpy import pi
from astropy import log
import os
try:
    from astropy.utils.console import ProgressBar
except ImportError:
    ProgressBar = lambda x: None
    ProgressBar.update = lambda x: None
import string
import struct

import time

def print_timing(func):
    """
    Prints execution time of decorated function.
    Included here because CLASS files can take a little while to read;
    this should probably be replaced with a progressbar
    """
    def wrapper(*arg,**kwargs):
        t1 = time.time()
        res = func(*arg,**kwargs)
        t2 = time.time()
        log.info('%s took %0.5g s' % (func.func_name, (t2-t1)))
        return res
    wrapper.__doc__ = func.__doc__
    return wrapper


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


There is also a nearly-complete specification in clic_file.f90


Empirical work on APEX mapping data:
    Index #383 is the last SgrA entry
    last index starts at 50048
    indexes seem to start at 1024
    each index is 128 bits
    there are 375 SGRA entries
    after the last SGRA entry, the next RAFGL entry is 252 bytes ahead

"""

""" Specification: http://iram.fr/IRAMFR/GILDAS/doc/html/class-html/node58.html """
filetype_dict = {'1A  ':'Multiple_IEEE','1   ':'Multiple_Vax','1B  ':'Multiple_EEEI',
    '9A  ':'Single_IEEE','9   ':'Single_Vax','9B  ':'Single_EEEI'}

record_lengths = {'1A': 512,
                  '2A': 1024*4}

header_id_numbers = {-2: 'GENERAL',
                     -3: 'POSITION',
                     -4: 'SPECTRO',
                     -5: 'BASELINE',
                    # -6: 'HISTORY',
                    # -8: 'SWITCH',
                     -10: 'DRIFT',
                     -14: 'CALIBRATION',
                    }

header_id_lengths = {-2: 9, # may really be 10?
                     -3: 17,
                     -4: 17,
                     -5: None, # variable length
                     -6: None, # variable length
                     -14: 25,
                    }


"""
GENERAL
 integer(kind=obsnum_length) :: num      ! [         ] Observation number
 integer(kind=4)             :: ver      ! [         ] Version number
 integer(kind=4)             :: teles(3) ! [         ] Telescope name
 integer(kind=4)             :: dobs     ! [MJD-60549] Date of observation
 integer(kind=4)             :: dred     ! [MJD-60549] Date of reduction
 integer(kind=4)             :: typec    ! [     code] Type of coordinates
 integer(kind=4)             :: kind     ! [     code] Type of data
 integer(kind=4)             :: qual     ! [     code] Quality of data
 integer(kind=4)             :: subscan  ! [         ] Subscan number
 integer(kind=obsnum_length) :: scan     ! [         ] Scan number
 ! Written in the entry
 real(kind=8)                :: ut       ! 1-2 [  rad] UT of observation
 real(kind=8)                :: st       ! 3-4 [  rad] LST of observation
 real(kind=4)                :: az       ! 5   [  rad] Azimuth
 real(kind=4)                :: el       ! 6   [  rad] Elevation
 real(kind=4)                :: tau      ! 7   [neper] Opacity
 real(kind=4)                :: tsys     ! 8   [    K] System temperature
 real(kind=4)                :: time     ! 9   [    s] Integration time
 ! Not in this section in file
 integer(kind=4)             :: xunit    ! [ code] X unit (if X coordinates section is present)
 ! NOT in data ---
 character(len=12)           :: cdobs    ! [string] Duplicate of dobs
 character(len=12)           :: cdred    ! [string] Duplicate of dred

"""

keys_lengths = {
        'unknown': [
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
     ],
        
        'GENERAL': [ # -2
     ( 'UT'      ,2,'float64'), #  rad UT of observation
     ( 'ST'      ,2,'float64'), #  rad LST of observation
     ( 'AZ'      ,1,'float32'), #  rad Azimuth
     ( 'EL'      ,1,'float32'), #  rad Elevation
     ( 'TAU'     ,1,'float32'), # neper Opacity
     ( 'TSYS'    ,1,'float32'), #    K System temperature
     ( 'TIME'    ,1,'float32'), #    s Integration time
                    # XUNIT should not be there?
     #( 'XUNIT'   ,1,'int32'),   # code X unit (if xcoord_sec is present)
     ] ,
     'POSITION': [ # -3
    ('SOURC',3,'|S12')  , #  [ ] Source name
    ('EPOCH',1,'float32'), #  [ ] Epoch of coordinates
    ('LAM'  ,2,'float64'), #[rad] Lambda
    ('BET'  ,2,'float64'), #[rad] Beta
    ('LAMOF',1,'float32'), #  [rad] Offset in Lambda
    ('BETOF',1,'float32'), #  [rad] Offset in Beta
    ('PROJ' ,1,'int32')  , # [rad] Projection system
    ('SL0P' ,1,'float64'), # lambda of descriptive system # MAY NOT EXIST IN OLD CLASS
    ('SB0P' ,1,'float64'), # beta of descriptive system   # MAY NOT EXIST IN OLD CLASS
    ('SK0P' ,1,'float64'), # angle of descriptive system  # MAY NOT EXIST IN OLD CLASS
    ],
     'SPECTRO': [ # -4
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
     'CALIBRATION': [ # -14
     ('ALIGN',1,'int32'),    # BUFFER (it's a zero - it is not declared in the docs!!!!)
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
     ('COUNT',3,'3float32'),   # [count] Power of Atm., Chopp., Cold
     ('LCALOF',1,'float32'),   # [ rad] Longitude offset for sky measurement
     ('BCALOF',1,'float32'),   # [ rad] Latitude offset for sky measurement
     ('GEOLONG',1,'float64'),   # [ rad] Geographic longitude of observatory # MAY NOT EXIST IN OLD CLASS
     ('GEOLAT',1,'float64'),   # [ rad] Geographic latitude of observatory   # MAY NOT EXIST IN OLD CLASS
         ],
    'BASELINE':[
        ('DEG',1,'int32'),       #! [         ] Degree of last baseline
        ('SIGFI',1,'float32'),   #! [Int. unit] Sigma
        ('AIRE',1,'float32'),    #! [Int. unit] Area under windows
        ('NWIND',1,'int32'),     #! [ ] Number of line windows
        # WARNING: These should probably have 'n', the second digit, = NWIND
        # The docs are really unclear about this, they say "W1(MWIND)"
        ('W1MWIND',1,'float32'), #! [km/s] Lower limits of windows
        ('W2MWIND',1,'float32'), #! [km/s] Upper limits of windows
        ('SINUS',3,'float32'),   #![]  Sinus baseline results
    ],

    'DRIFT':[ # 16?
        ('FREQ',1,'float64') ,  #! [ MHz] Rest frequency                   real(kind=8)    :: 
        ('WIDTH',1,'float32'),  #! [ MHz] Bandwidth                        real(kind=4)    :: 
        ('NPOIN',1,'int32')  ,  #! [    ] Number of data points              integer(kind=4) :: 
        ('RPOIN',1,'float32'),  #! [    ] Reference point                  real(kind=4)    :: 
        ('TREF',1,'float32') ,  #! [   ?] Time at reference                real(kind=4)    :: 
        ('AREF',1,'float32') ,  #! [ rad] Angular offset at ref.           real(kind=4)    :: 
        ('APOS',1,'float32') ,  #! [ rad] Position angle of drift          real(kind=4)    :: 
        ('TRES',1,'float32') ,  #! [   ?] Time resolution                  real(kind=4)    :: 
        ('ARES',1,'float32') ,  #! [ rad] Angular resolution               real(kind=4)    :: 
        ('BAD',1,'float32')  ,  #! [    ] Blanking value                   real(kind=4)    :: 
        ('CTYPE',1,'int32')  ,  #! [code] Type of offsets                    integer(kind=4) :: 
        ('CIMAG',1,'float64'),  #! [ MHz] Image frequency                  real(kind=8)    :: 
        ('COLLA',1,'float32'),  #! [   ?] Collimation error Az             real(kind=4)    :: 
        ('COLLE',1,'float32'),  #! [   ?] Collimation error El             real(kind=4)    :: 
    ],
     
     }

def _read_bytes(f, n):
    '''Read the next `n` bytes (from idlsave)'''
    return f.read(n)

"""
Warning: UNCLEAR what endianness should be!
Numpy seemed to get it right, and I think numpy assumes NATIVE endianness
"""

def _read_byte(f):
    '''Read a single byte (from idlsave)'''
    return numpy.uint8(struct.unpack('=B', f.read(4)[:1])[0])

def _read_int16(f):
    '''Read a signed 16-bit integer (from idlsave)'''
    return numpy.int16(struct.unpack('=h', f.read(4)[2:4])[0])

def _read_int32(f):
    '''Read a signed 32-bit integer (from idlsave)'''
    return numpy.int32(struct.unpack('=i', f.read(4))[0])

def _read_float32(f):
    '''Read a 32-bit float (from idlsave)'''
    return numpy.float32(struct.unpack('=f', f.read(4))[0])

def _align_32(f):
    '''Align to the next 32-bit position in a file (from idlsave)'''

    pos = f.tell()
    if pos % 4 != 0:
        f.seek(pos + 4 - pos % 4)
    return

def _read_word(f,length):
    if length > 0:
        chars = _read_bytes(f, length)
        _align_32(f)
    else:
        chars = None
    return chars

def _read_int(f):
    return struct.unpack('i',f.read(4))

def is_ascii(s):
    try:
        s.decode('ascii')
        return True
    except UnicodeDecodeError:
        return False
    except UnicodeEncodeError:
        return False


"""
from clic_file.f90: v1, v2
    integer(kind=4)  :: bloc       !  1   : observation address [records]       integer(kind=8)  :: bloc       !  1- 2: observation address [records]     integer(kind=4)   :: bloc     !  1   : block read from index
    integer(kind=4)  :: num        !  2   : observation number                  integer(kind=4)  :: word       !  3   : address offset      [4-bytes]     integer(kind=4)   :: num      !  2   : number read
    integer(kind=4)  :: ver        !  3   : observation version                 integer(kind=4)  :: ver        !  4   : observation version               integer(kind=4)   :: ver      !  3   : version read from index
    integer(kind=4)  :: sourc(3)   !  4- 6: source name                         integer(kind=8)  :: num        !  5- 6: observation number                character(len=12) :: csour    !  4- 6: source read from index
    integer(kind=4)  :: line(3)    !  7- 9: line name                           integer(kind=4)  :: sourc(3)   !  7- 9: source name                       character(len=12) :: cline    !  7- 9: line read from index
    integer(kind=4)  :: teles(3)   ! 10-12: telescope name                      integer(kind=4)  :: line(3)    ! 10-12: line name                         character(len=12) :: ctele    ! 10-12: telescope read from index
    integer(kind=4)  :: dobs       ! 13   : observation date    [class_date]    integer(kind=4)  :: teles(3)   ! 13-15: telescope name                    integer(kind=4)   :: dobs     ! 13   : date obs. read from index
    integer(kind=4)  :: dred       ! 14   : reduction date      [class_date]    integer(kind=4)  :: dobs       ! 16   : observation date    [class_date]  integer(kind=4)   :: dred     ! 14   : date red. read from index
    real(kind=4)     :: off1       ! 15   : lambda offset       [radian]        integer(kind=4)  :: dred       ! 17   : reduction date      [class_date]  real(kind=4)      :: off1     ! 15   : read offset 1
    real(kind=4)     :: off2       ! 16   : beta offset         [radian]        real(kind=4)     :: off1       ! 18   : lambda offset       [radian]      real(kind=4)      :: off2     ! 16   : read offset 2
    integer(kind=4)  :: typec      ! 17   : coordinates types                   real(kind=4)     :: off2       ! 19   : beta offset         [radian]      integer(kind=4)   :: type     ! 17   : type of read offsets
    integer(kind=4)  :: kind       ! 18   : data kind                           integer(kind=4)  :: typec      ! 20   : coordinates types                 integer(kind=4)   :: kind     ! 18   : type of observation
    integer(kind=4)  :: qual       ! 19   : data quality                        integer(kind=4)  :: kind       ! 21   : data kind                         integer(kind=4)   :: qual     ! 19   : Quality read from index
    integer(kind=4)  :: scan       ! 20   : scan number                         integer(kind=4)  :: qual       ! 22   : data quality                      integer(kind=4)   :: scan     ! 20   : Scan number read from index
    integer(kind=4)  :: proc       ! 21   : procedure type                      integer(kind=4)  :: scan       ! 23   : scan number                       real(kind=4)      :: posa     ! 21   : Position angle
    integer(kind=4)  :: itype      ! 22   : observation type                    integer(kind=4)  :: proc       ! 24   : procedure type                    integer(kind=4)   :: subscan  ! 22   : Subscan number
    real(kind=4)     :: houra      ! 23   : hour angle          [radian]        integer(kind=4)  :: itype      ! 25   : observation type                  integer(kind=4)   :: pad(10)  ! 23-32: Pad to 32 words
    integer(kind=4)  :: project    ! 24   : project name                        real(kind=4)     :: houra      ! 26   : hour angle          [radian]
    integer(kind=4)  :: pad1       ! 25   : unused word                         integer(kind=4)  :: project(2) ! 27   : project name
    integer(kind=4)  :: bpc        ! 26   : baseline bandpass cal status        integer(kind=4)  :: bpc        ! 29   : baseline bandpass cal status
    integer(kind=4)  :: ic         ! 27   : instrumental cal status             integer(kind=4)  :: ic         ! 30   : instrumental cal status
    integer(kind=4)  :: recei      ! 28   : receiver number                     integer(kind=4)  :: recei      ! 31   : receiver number
    real(kind=4)     :: ut         ! 29   : UT                  [s]             real(kind=4)     :: ut         ! 32   : UT                  [s] 
    integer(kind=4)  :: pad2(3)    ! 30-32: padding to 32 4-bytes word

equivalently

 integer(kind=obsnum_length) :: num      ! [         ] Observation number
 integer(kind=4)             :: ver      ! [         ] Version number
 integer(kind=4)             :: teles(3) ! [         ] Telescope name
 integer(kind=4)             :: dobs     ! [MJD-60549] Date of observation
 integer(kind=4)             :: dred     ! [MJD-60549] Date of reduction
 integer(kind=4)             :: typec    ! [     code] Type of coordinates
 integer(kind=4)             :: kind     ! [     code] Type of data
 integer(kind=4)             :: qual     ! [     code] Quality of data
 integer(kind=4)             :: subscan  ! [         ] Subscan number
 integer(kind=obsnum_length) :: scan     ! [         ] Scan number
"""

def _read_index(f, DEBUG=False, clic=False):
    x0 = f.tell()
    index = {
                "XBLOC":_read_int32(f),
                "XNUM":_read_int32(f),
                "XVER":_read_int32(f),
                "XSOURC":_read_word(f,12),
                "XLINE":_read_word(f,12),
                "XTEL":_read_word(f,12),
                "XDOBS":_read_int32(f),
                "XDRED":_read_int32(f),
                "XOFF1":_read_float32(f),# 		 first offset (real, radians) 
                "XOFF2":_read_float32(f),# 		 second offset (real, radians) 
                "XTYPE":_read_int32(f),# 		 coordinate system ('EQ'', 'GA', 'HO') 
                "XKIND":_read_int32(f),# 		 Kind of observation (0: spectral, 1: continuum, ) 
                "XQUAL":_read_int32(f),# 		 Quality (0-9)  
                "XSCAN":_read_int32(f),# 		 Scan number 
            }
    if clic: # use header set up in clic
        nextchunk = {
                    "XPROC":_read_int32(f),# "procedure type"
                    "XITYPE":_read_int32(f),#
                    "XHOURANG":_read_float32(f),#
                    "XPROJNAME":_read_int32(f),#
                    "XPAD1":_read_int32(f),
                    "XBPC" :_read_int32(f),
                    "XIC" :_read_int32(f),
                    "XRECEI" :_read_int32(f),
                    "XUT":_read_float32(f),
                    "XPAD2":numpy.fromfile(f,count=3,dtype='int32') # BLANK is NOT ALLOWED!!! It is a special KW
        }
    else:
        nextchunk = {"XPOSA":_read_float32(f),
                     "XSUBSCAN":_read_int32(f),
                     'XPAD2': numpy.fromfile(f,count=10,dtype='int32'),
                     }

    index.update(nextchunk)

    if f.tell() - x0 != 128:
        X = f.read(128-(f.tell()-x0))
        if DEBUG: print "read_index missed %i bits: %s" % (128-(f.tell()-x0),X)
        #raise IndexError("read_index did not successfully read 128 bytes at %i.  Read %i bytes." % (x0,f.tell()-x0))
    if any(not is_ascii(index[x]) for x in ('XSOURC','XLINE','XTEL')):
        raise ValueError("Invalid index read from {0}.".format(x0))
    return index

def _read_header(f,type=0):
    """
    Read a header entry from a CLASS file
    (helper function)
    """
    if type in keys_lengths:
        hdrsec = [(x[0],numpy.fromfile(f,count=1,dtype=x[2])[0])
                for x in keys_lengths[type]]
        return dict(hdrsec)
    pass


def _read_obshead(f,verbose=False):
    """
    Read the observation header of a CLASS file
    (helper function for read_class; should not be used independently)
    """
    #IDcode = f.read(4)
    #if verbose: print "IDcode: ",IDcode
    #if IDcode.strip() != '2':
    #    raise IndexError("Failure at %i" % (f.tell() - 4))
    nblocks,nbyteob,data_address,nheaders,data_length,obindex,nsec,obsnum= numpy.fromfile(f,count=8,dtype='int32')
    if verbose:
        print "nblocks,nbyteob,data_address,data_length,nheaders,obindex,nsec,obsnum",nblocks,nbyteob,data_address,data_length,nheaders,obindex,nsec,obsnum
        print "DATA_LENGTH: ",data_length
    
    seccodes = numpy.fromfile(f,count=nsec,dtype='int32')
    # Documentation says addresses then length: It is apparently wrong
    seclen = numpy.fromfile(f,count=nsec,dtype='int32')
    secaddr = numpy.fromfile(f,count=nsec,dtype='int32')
    if verbose: print "Section codes, addresses, lengths: ",seccodes,secaddr,seclen

    hdr = {'NBLOCKS':nblocks, 'NBYTEOB':nbyteob, 'DATAADDR':data_address,
            'DATALEN':data_length, 'NHEADERS':nheaders, 'OBINDEX':obindex,
            'NSEC':nsec, 'OBSNUM':obsnum}

    #return obsnum,seccodes
    return obsnum,hdr,dict(zip(seccodes,secaddr))

# THIS IS IN READ_OBSHEAD!!!
# def _read_preheader(f):
#     """
#     Not entirely clear what this is, but it is stuff that precedes the actual data
# 
#     Looks something like this:
#     array([          1,          -2,          -3,          -4,         -14,
#                  9,          17,          18,          25,          55,
#                 64,          81,          99, -1179344801,   979657591,
# 
#     -2, -3, -4, -14 indicate the 4 header types
#     9,17,18,25 *MAY* indicate the number of bytes in each
#     
# 
#     HOW is it indicated how many entries there are?
#     """
#     # 13 comes from counting 1, -2,....99 above
#     numbers = np.fromfile(f, count=13, dtype='int32')
#     sections = [n for n in numbers if n in header_id_numbers]
#     return sections

def downsample_1d(myarr,factor,estimator=np.mean):
    """
    Downsample a 1D array by averaging over *factor* pixels.
    Crops right side if the shape is not a multiple of factor.

    This code is pure numpy and should be fast.

    keywords:
        estimator - default to mean.  You can downsample by summing or
            something else if you want a different estimator
            (e.g., downsampling error: you want to sum & divide by sqrt(n))
    """
    if myarr.ndim != 1:
        raise ValueError("Only works on 1d data.  Says so in the title.")
    xs = myarr.size
    crarr = myarr[:xs-(xs % int(factor))]
    dsarr = estimator(np.concatenate([[crarr[i::factor] for i in
                                       range(factor)]]),axis=0)
    return dsarr


@print_timing
def read_class(filename,  DEBUG=False, apex=False, skip_blank_spectra=False,
               stop_position=None, skip_data=False,
               downsample_factor=None, memmap=False):
    """
    A hacked-together method to read a binary CLASS file.  It is strongly dependent on the incomplete
    `GILDAS CLASS file type Specification <http://iram.fr/IRAMFR/GILDAS/doc/html/class-html/node58.html>`_

    Parameters
    ----------
    filename: str
    skip_data: bool
        Read the file but don't store any of the data in memory.  Useful for
        getting headers.
    downsample_factor: None or int
        Factor by which to downsample data by averaging.  Useful for
        overresolved data.
    """
    t0 = time.time()
    statinfo = os.stat(filename)
    filelen = statinfo.st_size
    log.info("Loading file {0} with length {1} ({2} GB):".format(filename,filelen,filelen/(1024L**3)))
    f = open(filename,'rb')
    f.seek(0)
    if memmap:
        my_memmap = numpy.memmap(filename, offset=0, dtype='float32', mode='r')


    filetype = f.read(4)
    if filetype in filetype_dict:
        log.info("File %s is type %s" % (filename,filetype_dict[filetype]))
        record_length = record_lengths[filetype.strip()]
        version = 'v1' if '1' in filetype else 'v2'
    else:
        raise TypeError("File type error: %s." % filetype)
    
    nextblock,nindex,nex,nrecords = numpy.fromfile(f,count=4,dtype='int32')
    if DEBUG: print "nextblock,nindex,nex,nrecords",nextblock,nindex,nex,nrecords
    firstblock=numpy.fromfile(f,count=nex,dtype='int32') #[_read_byte(f) for ii in xrange(nex)]
    # purely empirical: HHT was always 512, apparently apex is 1024?
    # (though I had apex=True and it sort of worked maybe before?!)
    # 9/27/2014: figured out that yeah, this really should be the first block
    # (v1 files have 2 records for the 'header'; I think v2 do too?  Or their records are just longer?)
    f.seek(1024)

    if DEBUG: print "firstblock",firstblock
    if DEBUG: print "Done with header stuff at position %i" % f.tell()

    indexes = []
    #for ii in xrange(nindex):
    #    f.seek(128*(ii)+(firstblock[0]-1)*512)
    #    index = _read_index(f)
    #    # OLD DEBUG if index['XLINE'] not in ('HCOP(3-2)   ','N2HP(3-2)   '):
    #    # OLD DEBUG     raise Exception("Stopped at %i" % ii)
    #    indexes.append(index)

    #if f.tell() % 128 != 0:
    #    f.seek((f.tell()/128 + 1)*128)

    spectra = []
    header_list  = []
    spcount = 0
    jj = -1
    log.info("Pre-reading took {0} seconds".format(time.time()-t0))
    if not DEBUG: pb = ProgressBar(filelen)
    startpos = -1 # debug tool: make sure we are not in an infinite loop
    while f.tell() < filelen:
        jj += 1
        if not DEBUG: pb.update(f.tell())
        else: log.info("Iteration {0}".format(jj))
        if f.tell() == startpos:
            raise ValueError("Infinite Loop")
        startpos = f.tell()
        if stop_position is not None and f.tell() > stop_position:
            break
        IDcode = f.read(4)
        if IDcode == '\x00\x00\x00\x00':
            """ Skip over all blanks """
            f.seek(startpos)
            x = numpy.fromfile(f,count=128/4,dtype='int32')
            skipcount = 0
            # must allow for zero skipping!!
            pos = f.tell()
            while (x==0).all() and len(x)!=0:
                skipcount += 1
                pos = f.tell()
                x = numpy.fromfile(f,count=128/4,dtype='int32')
            f.seek(pos)
            if pos <= startpos:
                raise ValueError("Went backwards or did not advance.")
            if DEBUG: print "Loop %i: Skipped %i entries starting at %i.  pos=%i" % (jj, skipcount, startpos, pos)
            if jj >= filelen/128:
                # pdb here instead of excepting so that we can continue & repeat
                import pdb; pdb.set_trace()
                #raise Exception("Infinite Loop")
            continue
        elif IDcode.strip() != '2':
            f.seek(startpos)
            if DEBUG: print "IDCODE = ",IDcode," at pos ",startpos
            index = _read_index(f)
            if index['XLINE']: #numpy.any([L in index['XLINE'] for L in line]):
                if DEBUG: print "Found an index at %i: " % startpos,index
                indexes.append(index)
                continue
            else:
                print "Failure at %i: %i (/512=%i)" % (jj,f.tell(),f.tell()/512)
                raise Exception("Failure at %i: %i" % (jj,f.tell()))
        else:
            if DEBUG: print "Reading HEADER and OBSERVATION at %i" % f.tell()
            obsnum,obshead,sections = _read_obshead(f,verbose=DEBUG)
        pos = f.tell()
        # first one of these starts
        if DEBUG:
            # Actually, it's pretty clear that this is leading in to the start of headers.
            print " UNCLEAR INFO AT %i PAST %i: %i" % (f.tell()-startpos,startpos,pos)
            junk =  f.read(168) 
            print "Junk as str: ",junk
            print "Junk as int: ",[struct.unpack('=i',junk[i*4:i*4+4])[0] for i in range(len(junk)/4)]
            print "Junk as flt: ",[struct.unpack('=f',junk[i*4:i*4+4])[0] for i in range(len(junk)/4)]
            #f.seek(startpos)
            #print numpy.fromfile(f,count=168/4+10,dtype='int32') 
            #f.seek(startpos)
            #print numpy.fromfile(f,count=168/4+10,dtype='float32') 
            #f.seek(startpos)
            #print numpy.fromfile(f,count=168+40,dtype='int8') 

        if apex: # there are 33 bytes of junk...
            # This is only for old, not modern, APEX files (and I don't know where the cutoff is)
            somejunk = f.read(33*4)
            print "APEX Junk as str: ",somejunk
            print "APEX Junk as int: ",[struct.unpack('=i',somejunk[i*4:i*4+4])[0] for i in range(len(somejunk)/4)]
            print "APEX Junk as flt: ",[struct.unpack('=f',somejunk[i*4:i*4+4])[0] for i in range(len(somejunk)/4)]
            f.seek(pos+84+12*4)
        else: # ????
            if filetype.strip() == '9A':
                if f.tell() != pos + 156:
                    #print "Wrong position %i, skipping to %i" % (f.tell(),pos+168)
                    f.seek(pos+156)
                    #warnings.warn("Skip from %i to %i" % (pos,pos+168))
            else:
                if f.tell() != pos + 168:
                    #print "Wrong position %i, skipping to %i" % (f.tell(),pos+168)
                    f.seek(pos+168)
                    #warnings.warn("Skip from %i to %i" % (pos,pos+168))

        header = obshead


        # datastart must be at the end of all these sections
        datastart = 0

        if DEBUG:
            print "Starting header reading at startpos=%i" % startpos

        for section_id,section_address in sections.iteritems():
            # Section addresses are 1-indexed byte addresses
            # in the current "block"
            f.seek(startpos + (section_address-1)*4)
            temp_hdr = _read_header(f, type=header_id_numbers[section_id])
            header.update(temp_hdr)
            datastart = max(datastart,f.tell())

        if DEBUG > 1:
            raise ValueError("Debug Breakpoint")

        # can't guarantee that the loop above goes in the right order,
        # so make sure we end up at the right spot
        f.seek(datastart)

        #Header2 = _read_header(f,type='GENERAL')
        #Header3 = _read_header(f,type='POSITION')
        #Header4 = _read_header(f,type='SPECTRO')
        #if DEBUG: print "Line %i (byte %i) - OBSERVATION %i (%i): %s, %s" % (f.tell(),(f.tell()-startpos)/4,obsnum,spcount,Header3['SOURC'],Header4['LINE']),
        #Header14 = _read_header(f,type='CALIBRATION')
        #if DEBUG: print "\nLine %i (byte %i) - CALIBRATION:" % (f.tell(),(f.tell()-startpos)/4)

        if not all((a in string.printable for a in header['SOURC'])):
            raise ValueError("Source has weird name.")

        ## purely empirical; NO idea why
        #if apex:
        #    f.seek(238*4+f.tell())

        #hdr = Header3
        #hdr.update(Header2)
        #hdr.update(Header4)
        #hdr.update(Header14)
        #hdr.update(unclear)
        hdr = header
        hdr.update(obshead) # re-overwrite things
        hdr.update({'OBSNUM':obsnum,'RECNUM':spcount})
        hdr.update({'RA':hdr['LAM']/pi*180,'DEC':hdr['BET']/pi*180})
        hdr.update({'RAoff':hdr['LAMOF']/pi*180,'DECoff':hdr['BETOF']/pi*180})
        hdr.update({'OBJECT':hdr['SOURC'].strip()})
        hdr.update({'BUNIT':'Tastar'})
        hdr.update({'EXPOSURE':hdr['TIME']})
        #for key in hdr:
        #    if pyfits.Card._comment_FSC_RE.match(str(hdr[key])) is None:
        #        print "Setting hdr[%s] to ''" % (key)
        #        hdr[key] = ''
        spcount += 1

        if -10 in sections and 'NPOIN' in hdr:
            npoin = hdr['NPOIN']
            nchan = npoin
        else:
            nchan = hdr['NCHAN']

        if nchan != hdr['DATALEN'] and not skip_blank_spectra:
            print "WARNING: nchan=%i, datalen=%i" % (nchan, hdr['NCHAN'])
            #raise ValueError("data_length != nchan")
        elif nchan != hdr['DATALEN'] and skip_blank_spectra:
            print "Loop %i: Skipped a spectrum with nchan=%i, datalen=%i at %i" % (jj, nchan, hdr['DATALEN'], f.tell())
            continue
        else:
            if DEBUG:
                print "Everything is fine.  %i" % jj
        if DEBUG:
            print "Spectrum has %i channels at %i" % (nchan,f.tell())
        if skip_data:
            # Skip the data: 4 bytes per float32
            f.seek(f.tell() + nchan*4)
        else:
            if memmap:
                here = f.tell()
                #spectrum = numpy.memmap(filename, offset=here, dtype='float32',
                #                        mode='r', shape=(nchan,))
                spectrum = my_memmap[here/4:here/4+nchan]
                f.seek(here+nchan*4)
            else:
                spectrum = numpy.fromfile(f,count=nchan,dtype='float32')
            if DEBUG > 2:
                print "First digits of spectrum: ",spectrum[:10]
            if downsample_factor is not None:
                spectra.append(downsample_1d(spectrum, downsample_factor))
            else:
                spectra.append(spectrum)

        if downsample_factor is not None:
            if 'NCHAN' in hdr:
                hdr['NCHAN'] /= downsample_factor
            elif 'NPOIN' in hdr:
                hdr['NPOIN'] /= downsample_factor
            else:
                log.info("Did not find any header keywords reporting NCHAN")
            if 'DATALEN' in hdr:
                hdr['DATALEN'] /= downsample_factor
            for kw in ['FRES','VRES']:
                if kw in hdr:
                    hdr[kw] *= downsample_factor
            if 'RCHAN' in hdr:
                scalefactor = 1./downsample_factor
                hdr['RCHAN'] = (hdr['RCHAN']-1)*scalefactor + 0.5 + scalefactor/2.
        header_list.append(hdr)
        if f.tell() % 128 != 0:
            discard = f.read(128 - f.tell() % 128)
            if DEBUG > 3:
                print "DISCARD: ",discard," POSITION: ",f.tell()
            #raise ValueError("Bad position : %i "%f.tell())
        #f.seek((f.tell()/nchan + 1)*nchan)
        if DEBUG:
            print "Continuing to sp# ",spcount+1," starting at ",f.tell()," out of ",filelen

    f.close()
    return spectra,header_list,indexes

def make_axis(header,imagfreq=False):
    """
    Create a :class:`pyspeckit.spectrum.units.SpectroscopicAxis` from the CLASS "header"
    """
    from .. import units

    rest_frequency = header.get('RESTF')
    xunits = 'MHz'
    nchan = header.get('NCHAN')
    voff = header.get('VOFF')
    foff = header.get('FOFF')
    doppler = header.get('DOPPLER')
    fres = header.get('FRES')
    refchan = header.get('RCHAN')
    imfreq = header.get('IMAGE')

    if not imagfreq:
        xarr =  rest_frequency + (numpy.arange(1, nchan+1) - refchan) * fres
        XAxis = units.SpectroscopicAxis(xarr,'MHz',frame='rest',refX=rest_frequency)
    else:
        xarr = imfreq - (numpy.arange(1, nchan+1) - refchan) * fres
        XAxis = units.SpectroscopicAxis(xarr,'MHz',frame='rest',refX=imfreq)

    return XAxis
    
import pyspeckit
@print_timing
def class_to_obsblocks(filename, telescope, line, datatuple=None, source=None,
                       imagfreq=False, DEBUG=False,  **kwargs):
    """
    Load an entire CLASS observing session into a list of ObsBlocks based on
    matches to the 'telescope', 'line' and 'source' names

    Parameters
    ----------
    filename : string
        The Gildas CLASS data file to read the spectra from.
    telescope : list
        List of telescope names to be matched.
    line : list
        List of line names to be matched.
    source : list (optional)
        List of source names to be matched. Defaults to None.
    imagfreq : bool
        Create a SpectroscopicAxis with the image frequency.
    """
    if datatuple is None:
        spectra,header,indexes = read_class(filename,DEBUG=DEBUG, **kwargs)
    else:
        spectra,header,indexes = datatuple

    obslist = []
    lastscannum = -1
    spectrumlist = None
    for sp,hdr,ind in zip(spectra,header,indexes):
        hdr.update(ind)
        # this is slow but necessary...
        H = pyfits.Header()
        for k,v in hdr.iteritems():
            if hasattr(v,"__len__") and not isinstance(v,str):
                if len(v) > 1:
                    for ii,vv in enumerate(v):
                        H.update(k[:7]+str(ii),vv)
                else:
                    H.update(k,v[0])
            elif pyfits.Card._comment_FSC_RE.match(str(v)) is not None:
                H.update(k,v)
        scannum = hdr['XSCAN']
        if 'XTEL' in hdr and hdr['XTEL'].strip() not in telescope:
            continue
        if hdr['LINE'].strip() not in line:
            continue
        if (source is not None) and (hdr['SOURC'].strip() not in source):
            continue
        hdr.update({'RESTFREQ':hdr.get('RESTF')})
        H.update('RESTFREQ',hdr.get('RESTF'))

        #print "Did not skip %s,%s.  Scannum, last: %i,%i" % (hdr['XTEL'],hdr['LINE'],scannum,lastscannum)

        if scannum != lastscannum:
            lastscannum = scannum
            if spectrumlist is not None:
                obslist.append(pyspeckit.ObsBlock(spectrumlist))
            xarr = make_axis(hdr,imagfreq=imagfreq)
            spectrumlist = [(
                pyspeckit.Spectrum(xarr=xarr,
                    header=H,
                    data=sp))]
        else:
            spectrumlist.append(
                pyspeckit.Spectrum(xarr=xarr,
                    header=H,
                    data=sp))

    return obslist

@print_timing
def class_to_spectra(filename, datatuple=None, **kwargs):
    """
    Load each individual spectrum within a CLASS file into a list of Spectrum
    objects
    """
    if datatuple is None:
        spectra,header,indexes = read_class(filename, **kwargs)
    else:
        spectra,header,indexes = datatuple

    spectrumlist = []
    for sp,hdr,ind in zip(spectra,header,indexes):
        hdr.update(ind)
        xarr = make_axis(hdr)
        spectrumlist.append(
            pyspeckit.Spectrum(xarr=xarr,
                               header=hdr,
                               data=sp))

    return pyspeckit.Spectra(spectrumlist)

def tests():
    """
    Tests are specific to the machine on which this code was developed. 
    """
    fn1 = '/Users/adam/work/bolocam/hht/class_003.smt' 
    #fn1 = '/Users/adam/work/bolocam/hht/class_001.smt' 
    #fn1 = '/Users/adam/work/bolocam/hht/test_SMT-F1M-VU-20824-073.cls' 
    #fn2 = '/Users/adam/work/bolocam/hht/test_SMT-F1M-VU-79472+203.cls'
    #F1 = read_class(fn1)#,DEBUG=True)
    #F2 = read_class(fn2)
    n2hp = class_to_obsblocks(fn1,telescope=['SMT-F1M-HU','SMT-F1M-VU'],line=['N2HP(3-2)','N2H+(3-2)'])
    hcop = class_to_obsblocks(fn1,telescope=['SMT-F1M-HL','SMT-F1M-VL'],line=['HCOP(3-2)','HCO+(3-2)'])
