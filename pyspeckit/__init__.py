"""
:Author: Adam Ginsburg <adam.g.ginsburg@gmail.com> and Jordan Mirocha <mirochaj@gmail.com>
.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>
"""
__all__ = ['spectrum','cubes','wrappers']
import spectrum
try:
    import cubes
except ImportError:
    print "pyspeckit.cubes module not imported - cubes requires pywcs"
import wrappers
from wrappers import *
from cubes import *
from spectrum import *
