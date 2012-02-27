"""
:Author: Adam Ginsburg <adam.g.ginsburg@gmail.com> and Jordan Mirocha <mirochaj@gmail.com>
.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>
"""
__all__ = ['spectrum','cubes','wrappers']
__version__ = '0.1.5'
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
