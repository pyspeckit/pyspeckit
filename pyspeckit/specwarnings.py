"""
Override components of python's builting warnings (as suggested by the manual)
"""

import warnings
import sys

class PyspeckitWarning(Warning):
    """
    The base warning class from which all pyspeckit warnings should inherit.
    """

def showwarning(message, category, filename, lineno, file=None, line=None):
    """Hook to write a warning to a file; replace if you like."""
    if file is None:
        file = sys.stderr
    try:
        msg = str(message) # tried using message.message before, but led to infinite recursion for deprecation warning
        file.write(msg+"\n")
        #file.write(formatwarning(message, category, filename, lineno, line))
    except IOError:
        pass # the file (probably stderr) is invalid - this warning gets lost.

# This is not overwritten, but may be in the future... kept here so we know what to edit
#def formatwarning(message, category, filename, lineno, line=None):
#    """Function to format a warning the standard way."""
#    s =  "%s:%s: %s: %s\n" % (filename, lineno, category.__name__, message)
#    line = linecache.getline(filename, lineno) if line is None else line
#    if line:
#        line = line.strip()
#        s += "  %s\n" % line
#    return s

warnings.showwarning = showwarning

warn = warnings.warn
