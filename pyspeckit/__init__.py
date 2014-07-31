"""
:Author: Adam Ginsburg <adam.g.ginsburg@gmail.com> and Jordan Mirocha <mirochaj@gmail.com>
.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>
"""
__all__ = ['spectrum','cubes','wrappers']
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
    execfile("../tests/run_tests.py")
