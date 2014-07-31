"""
:Author: Adam Ginsburg <adam.g.ginsburg@gmail.com> and Jordan Mirocha <mirochaj@gmail.com>
.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>
"""
from ._astropy_init import *

__all__ = ['spectrum','cubes','wrappers']
if not _ASTROPY_SETUP_:
    from version import version as __version__
    import spectrum
    import specwarnings
    try:
        import cubes
    except ImportError:
        specwarnings.warn( "pyspeckit.cubes module not imported - cubes requires pywcs" )
    import wrappers
    from wrappers import *
    from cubes import *
    from spectrum import *

    def test(*args, **kwargs):
        import tests
        #import os
        #os.chdir(os.path.split(os.path.abspath(tests.__file__))[0])
        tests.run_tests.test_everything()
