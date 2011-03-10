"""
writers.py

Author: Jordan Mirocha
Affiliation: University of Colorado at Boulder
Created on 2011-03-09.

Description: Write out fitting results to ascii or hdf5.     
"""

import os, os.path
from config import *

h5fail = False
try: import h5py
except ImportError: h5fail = True

#if 'hdf5' in writer_dict: writer = HDF5Writer else: print "ERROR: no hdf5 writer installed"

class Writer(object):
    def __init__(self, Spectrum):
        self.Spectrum = Spectrum
        self.cfg = spcfg.cfg
        
        if spcfg.cfg['data_format'] in ['hdf5', 'h5'] and h5fail: 
            print 'Failed to import h5py, will write out in ascii format instead.'
            self.write_data = self.write_ascii
        elif spcfg.cfg['data_format'] in ['hdf5', 'h5']: 
            self.write_data = self.write_hdf5
        else: 
            self.write_data = self.write_ascii
        
    def __call__(self):
        self.write_data()    

    def write_hdf5(self, clobber = True):
        """
        Write all fit information to hdf5 file.
        """
        
        fn = "{0}_fit.h5".format(self.Spectrum.fileprefix)
        if not clobber:
            i = 1
            while os.path.exists(fn):
                fn = "{0}_fit({1}).h5".format(self.Spectrum.fileprefix, i)
                i += 1
                
        f = h5py.File(fn, 'w')
        
        f.create_dataset('xarr', data = self.Spectrum.xarr[self.Spectrum.specfit.gx1:self.Spectrum.specfit.gx2])
        f.create_dataset('model_spectrum', data = self.Spectrum.specfit.model)
        f.create_dataset('model_params', data = self.Spectrum.specfit.modelpars)
        f.create_dataset('model_errs', data = self.Spectrum.specfit.modelerrs)
        f.create_dataset('residuals', data = self.Spectrum.specfit.residuals)
        for i, element in enumerate(self.Spectrum.specfit.modelcomponents):
            f.create_dataset('model_component{0}'.format(i), data = self.Spectrum.specfit.modelcomponents[i])
            
        f.close()
        
    def write_ascii(self, clobber = True):
        """
        Write all fit information to an ASCII file.
        """
        
        fn = "{0}_fit.dat".format(self.Spectrum.fileprefix)
        if not clobber:
            i = 1
            while os.path.exists(fn):
                fn = "{0}_fit({1}).dat".format(self.Spectrum.fileprefix, i)
                i += 1
                
        f = open(fn, 'w')
        
        # Print header
        print >> f, "# Column 1: {0}".format("x-values")
        print >> f, "# Column 2: {0}".format("model spectrum")
        for i, element in enumerate(self.Spectrum.specfit.modelcomponents):
            print >> f, "# Column {0}: model spectrum component {1}".format(i + 3, i + 1)        
        print >> f, "# Column {0}: residuals".format(i + 4)
        print >> f, ""    
                
        components = zip(*self.Spectrum.specfit.modelcomponents)        
        for i, element in enumerate(self.Spectrum.specfit.model):
            line = "{0:10}{1:10}".format(self.Spectrum.xarr[self.Spectrum.specfit.gx1:self.Spectrum.specfit.gx2][i], 
                round(self.Spectrum.specfit.model[i], 5))
            for j, component in enumerate(components[i]): line += "{0:10}".format(round(component, 5))    
            line += "{0:10}".format(round(self.Spectrum.specfit.residuals[i], 5))       
                
            print >> f, line
            
        print >> f, ""
            
        f.close()
        
        
        
        
        
        
        
        
        
