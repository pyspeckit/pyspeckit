writers = {}
suffix_types = {}

class Writer(object):
    def __init__(self, Spectrum):
        self.Spectrum = Spectrum
        
    def __call__(self,**kwargs):
        self.write_data(**kwargs)

    def write_data(self):
        print "This is a dummy function intended to be overloaded."

import hdf5_writer
write_hdf5 = hdf5_writer.write_hdf5
import fits_writer
write_fits = fits_writer.write_fits
