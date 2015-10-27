"""
"""
from .classes import Spectrum,Spectra,ObsBlock
from .OpticalSpectrum import OpticalSpectrum
from . import fitters,plotters,baseline,units
from . import smooth
from . import correlate
from . import headers
from . import moments
from . import units
from . import utils
from . import readers
from . import writers
from . import logger
from .. import config

def register_reader(filetype, function, suffix, default=False):
    ''' 
    Register a reader function.

    Parameters
    ----------
    filetype: str
        The file type name
    function: function
        The reader function.  Should take a filename as input and return an
        X-axis object (see units.py), a spectrum, an error spectrum (initialize
        it to 0's if empty), and a pyfits header instance
    suffix: int
        What suffix should the file have?


    '''

    readers.readers[filetype] = function
    if suffix in readers.suffix_types:
        if default:
            readers.suffix_types[suffix].insert(0,filetype)
        else:
            readers.suffix_types[suffix].append(filetype)
    else: # if it's the first, it defaults to default!
        readers.suffix_types[suffix] = [filetype]

register_reader('fits',readers.open_1d_fits,'fits',default=True)
register_reader('fits',readers.open_1d_fits,'fit')
register_reader('sdss',readers.read_sdss,'fits')
register_reader('pyfits',readers.open_1d_pyfits,'')
register_reader('txt',readers.open_1d_txt,'txt')
register_reader('txt',readers.open_1d_txt,'dat')
register_reader('tspec',readers.tspec_reader,'fits')
register_reader('hdf5',readers.open_hdf5,'hdf5')

def register_writer(filetype, function, suffix, default=False):
    ''' 
    Register a writer function.

    Parameters
    ----------
    filetype:  string 
        The file type name
    function:  function 
        The writer function.  Will be an attribute of Spectrum
        object, and called as spectrum.Spectrum.write_hdf5(), 
        for example. 
    suffix:  int 
        What suffix should the file have?

    '''

    writers.writers[filetype] = function
    if suffix in writers.suffix_types:
        if default:
            writers.suffix_types[suffix].insert(0,filetype)
        else:
            writers.suffix_types[suffix].append(filetype)
    else: # if it's the first, it defaults to default!
        writers.suffix_types[suffix] = [filetype]

register_writer('hdf5',writers.write_hdf5,'hdf5',default=False)
register_writer('fits',writers.write_fits,'fits',default=True)
#register_writer('txt',writers.write_ascii,'txt')
#register_writer('txt',writers.write_ascii,'txt')
