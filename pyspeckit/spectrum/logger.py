"""
Logger for the spectroscopic toolkit packge

Goal: Use decorators to log each command in full

Author: Adam Ginsburg
Created: 03/17/2011
"""
from __future__ import print_function
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
        print("Began logging at %s" % (time.strftime("%M/%d/%Y %H:%M:%S",time.localtime())), file=self.outfile)

    def __call__(self, function):
        """
        """
        def wrapper(*args, **kwargs):
            modname = function.__module__
            fname = function.__name__
            print("%s.%s" % (modname,fname) +\
                    "("+"".join([a.__name__ for a in args]) + \
                    "".join(["%s=%s" % (k,v) for (k,v) in kwargs])+ ")", file=self.outfile)
            return function
        return wrapper

    def close(self):
        self.outfile.close()
