"""
pyspeckit is a package for spectrum plotting, fitting, and analysis

For details, see https://pyspeckit.readthedocs.io/
"""

__all__ = ['spectrum','cubes','wrappers']

from . import spectrum
from . import specwarnings
from . import cubes
from . import wrappers
from .wrappers import *
from .cubes import *
from .spectrum import *

# Import version from version module
try:
    from .version import __version__
except ImportError:
    # Fallback version if version module is not available
    __version__ = '1.0.5.dev'

try:
    from .tests import run_tests

    def test(*args, **kwargs):
        #import os
        #os.chdir(os.path.split(os.path.abspath(tests.__file__))[0])
        from .spectrum.tests import test_eqw
        test_eqw.test_eqw()
        from .spectrum.models.tests import test_template
        test_template.test_template()
        test_template.test_template_withcont()
        run_tests.test_everything()
except ImportError:
    # This makes no sense.
    pass
