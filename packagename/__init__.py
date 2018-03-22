# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

# Uncomment to enforce Python version check during package import.
# This is the same check as the one at the top of setup.py
#class UnsupportedPythonError(Exception):
#    pass
#__minimum_python_version__ = '3.5'
#if sys.version_info < tuple((int(val) for val in __minimum_python_version__.split('.'))):
#    raise UnsupportedPythonError("{} does not support Python < {}".format(__package__, __minimum_python_version__))

if not _ASTROPY_SETUP_:
    # For egg_info test builds to pass, put package imports here.

    from .example_mod import *
