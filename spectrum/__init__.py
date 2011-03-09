import readers,fitters,plotters,baseline
from spectrum import Spectrum
import units

def register_fitter(name, function, npars, multisingle='single',
        override=False):
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
    '''

    if multisingle == 'single':
        if not name in fitters.singlefitters or override:
            fitters.singlefitters[name] = function
    elif multisingle == 'multi':
        if not name in fitters.multifitters or override:
            fitters.multifitters[name] = function
    else:
        raise Exception("Fitting function %s is already defined" % name)

    fitters.npars[name] = npars


import models
register_fitter('ammonia',models.ammonia_model(),6,multisingle='multi')
register_fitter('gaussian',models.gaussian_fitter(multisingle='multi'),3,multisingle='multi')
register_fitter('gaussian',models.gaussian_fitter(multisingle='single'),3,multisingle='single')
register_fitter('voigt',models.voigt_fitter(multisingle='multi'),4,multisingle='multi')
register_fitter('voigt',models.voigt_fitter(multisingle='single'),4,multisingle='single')
