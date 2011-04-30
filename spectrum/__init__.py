import fitters,plotters,baseline,units
from spectrum import Spectrum,Spectra,ObsBlock
import smooth
import logger
import config

def register_fitter(name, function, npars, multisingle='single',
        override=False, key=None):
    ''' 
    Register a fitter function.

    Required Arguments:

        *name*: [ string ]
            The fit function name. 

        *function*: [ function ]
            The fitter function.  Single-fitters should take npars + 1 input
            parameters, where the +1 is for a 0th order baseline fit.  They
            should accept an X-axis and data and standard fitting-function
            inputs (see, e.g., gaussfitter).  Multi-fitters should take N *
            npars, but should also operate on X-axis and data arguments.

        *npars*: [ int ]
            How many parameters does the function being fit accept?

    Optional Keyword Arguments:

        *multisingle*: [ 'multi' | 'single' ] 
            Is the function a single-function fitter (with a background), or
            does it allow N copies of the fitting function?

        *override*: [ True | False ]
            Whether to override any existing type if already present.

        *key*: [ char ]
            Key to select the fitter in interactive mode
    '''

    if multisingle == 'single':
        if not name in fitters.singlefitters or override:
            fitters.singlefitters[name] = function
    elif multisingle == 'multi':
        if not name in fitters.multifitters or override:
            fitters.multifitters[name] = function
    elif name in fitters.singlefitters or name in fitters.multifitters:
        raise Exception("Fitting function %s is already defined" % name)

    if key is not None:
        fitters.fitkeys[key] = name
        fitters.interactive_help_message += "\n'%s' - select fitter %s" % (key,name)
    fitters.npars[name] = npars


import models
register_fitter('ammonia',models.ammonia_model(),6,multisingle='multi',key='a')
register_fitter('gaussian',models.gaussian_fitter(multisingle='multi'),3,multisingle='multi',key='g')
register_fitter('gaussian',models.gaussian_fitter(multisingle='single'),3,multisingle='single')
register_fitter('voigt',models.voigt_fitter(multisingle='multi'),4,multisingle='multi',key='v')
register_fitter('voigt',models.voigt_fitter(multisingle='single'),4,multisingle='single')

import readers
def register_reader(filetype, function, suffix, default=False):
    ''' 
    Register a reader function.

    Required Arguments:

        *filetype*: [ string ]
            The file type name

        *function*: [ function ]
            The reader function.  Should take a filename as input
            and return an X-axis object (see units.py), a spectrum,
            an error spectrum (initialize it to 0's if empty), and a
            pyfits header instance

        *suffix*: [ int ]
            What suffix should the file have?

    Optional Keyword Arguments:

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
register_reader('pyfits',readers.open_1d_pyfits,'')
register_reader('txt',readers.open_1d_txt,'txt')
register_reader('txt',readers.open_1d_txt,'dat')
register_reader('tspec',readers.tspec_reader,'fits')
register_reader('hdf5',readers.open_hdf5,'hdf5')

import writers
def register_writer(filetype, function, suffix, default=False):
    ''' 
    Register a writer function.

    Required Arguments:

        *filetype*: [ string ]
            The file type name

        *function*: [ function ]
            The writer function.  Will be an attribute of Spectrum
            object, and called as spectrum.Spectrum.write_hdf5(), 
            for example. 

        *suffix*: [ int ]
            What suffix should the file have?

    Optional Keyword Arguments:

    '''

    writers.writers[filetype] = function
    if suffix in writers.suffix_types:
        if default:
            writers.suffix_types[suffix].insert(0,filetype)
        else:
            writers.suffix_types[suffix].append(filetype)
    else: # if it's the first, it defaults to default!
        writers.suffix_types[suffix] = [filetype]

# Currently, names of writers must be unique! (That should probably always be true)
register_writer('hdf5',writers.write_hdf5,'hdf5',default=False)
register_writer('fits',writers.write_fits,'fits',default=True)
#register_writer('txt',writers.write_ascii,'txt')
#register_writer('txt',writers.write_ascii,'txt')
