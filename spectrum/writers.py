"""
writers.py

Author: Jordan Mirocha
Affiliation: University of Colorado at Boulder
Created on 2011-08-09.

Description: Write out fitting results to ascii or hdf5.     
"""

import os, os.path
from config import *

h5fail = False
try: import h5py
except ImportError: h5fail = True

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
        
        
