"""
History logger for the spectroscopic toolkit packge

Goal: Save history...
(should this be implemented as a Header class method?)

Author: Adam Ginsburg
Created: 03/18/2011
"""
import time
import pyfits

def write_history(header, string):
    """
    Add a line to the header's history
    """
    if isinstance(header, pyfits.Header):
        hdrstring = time.strftime("%m/%d/%y %H:%M:%S",time.localtime()) + " " + string
        try:
            header.add_history(hdrstring)
        except AttributeError:
            print "WARNING: Error in history writing.  Could not add this string: %s" % hdrstring

    else:
        raise TypeError("Error: header is not a header instance!")
