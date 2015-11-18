from __future__ import print_function
"""
History logger for the spectroscopic toolkit packge

Goal: Save history...
(should this be implemented as a Header class method?)

Author: Adam Ginsburg
Created: 03/18/2011
"""
import time
from ..specwarnings import warn
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits

def write_history(header, string):
    """
    Add a line to the header's history
    """
    if not isinstance(header, pyfits.Header):
        warn("header is not a header instance!")
    hdrstring = time.strftime("%m/%d/%y %H:%M:%S",time.localtime()) + " " + string
    try:
        header.add_history(hdrstring)
    except AttributeError:
        print("WARNING: Error in history writing.  Could not add this string: %s" % hdrstring)
