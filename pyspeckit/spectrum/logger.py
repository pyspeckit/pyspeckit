"""
Logger for the spectroscopic toolkit packge

Goal: Use decorators to log each command in full

Author: Adam Ginsburg
Created: 03/17/2011
"""
import os
import time

class Logger(object):
    """
    Logger object.  Should be initiated on import.
    """

    def __init__(self, filename):
        """
        Open a log file for writing
        """

        newfilename = filename
        ii=0
        while os.path.exists(filename):
            newfilename = filename+".%3i" % ii
            ii += 1
        if newfilename != filename:
            os.rename(filename,newfilename)

        self.outfile = open(filename,'w')
        print >>self.outfile,"Began logging at %s" % (time.strftime("%M/%d/%Y %H:%M:%S",time.localtime()))

    def __call__(self, function):
        """
        """
        def wrapper(*args, **kwargs):
            modname = function.__module__
            fname = function.__name__
            print >>self.outfile,"%s.%s" % (modname,fname) +\
                    "("+"".join([a.__name__ for a in args]) + \
                    "".join(["%s=%s" % (k,v) for (k,v) in kwargs])+ ")"
            return function
        return wrapper

    def close(self):
        self.outfile.close()
