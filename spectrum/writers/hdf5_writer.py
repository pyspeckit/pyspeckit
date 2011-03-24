import os

h5check = True
try: import h5py
except ImportError: h5check = False

class write_hdf5(object):
    def __init__(self, Spectrum):
        self.Spectrum = Spectrum
        
    def write_data(self, filename = None, extras = None, clobber = True):
        """
        Write information to hdf5 file.
        """
        
        if not h5check: print "Cannot write to hdf5 - h5py import failed."
        
        else: 
            if filename is None: fn = "{0}_out.h5".format(self.Spectrum.fileprefix)
            else: fn = filename
            
            if clobber: f = h5py.File(fn, 'w')
            else:
                i = 1
                while os.path.exists(fn):
                    fn = "{0}_fit({1}).h5".format(self.Spectrum.fileprefix, i)
                    i += 1
                            
            # By default, write out x-array, spectrum, errors
            f.create_dataset('xarr', data = self.Spectrum.xarr)
            f.create_dataset('data', data = self.Spectrum.data)
            f.create_dataset('err', data = self.Spectrum.error)
            
            if extras is not None:
                for extra in extras:
                    f.create_dataset(extras[0], data = extras[1])
                
            f.close()