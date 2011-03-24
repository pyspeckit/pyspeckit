import os

h5check = True
try: import h5py
except ImportError: h5check = False

class write_hdf5(object):
    def __init__(self, Spectrum):
        self.Spectrum = Spectrum
        
    def write_data(self, filename = None, newsuffix = 'out', extras = None, clobber = True):
        """
        Write information to hdf5 file.
        
        extras = [(dataset name, data), (dataset name, data)] 
        
        To do: leave option of writing to groups (for model components especially?)
        """
        
        if not h5check: print "Cannot write to hdf5 - h5py import failed."
        
        else: 
            if filename is None: fn = "{0}_{1}.hdf5".format(self.Spectrum.fileprefix, newsuffix)
            else: fn = filename
            
            if clobber: f = h5py.File(fn, 'w')
            else:
                i = 1
                while os.path.exists(fn):
                    fn = "{0}_fit({1}).hdf5".format(self.Spectrum.fileprefix, i)
                    i += 1
                            
            # By default, write out x-array, spectrum, errors
            xarr = f.create_dataset('xarr', data = self.Spectrum.xarr)
            data = f.create_dataset('data', data = self.Spectrum.data)
            error = f.create_dataset('error', data = self.Spectrum.error)
            
            # Add metadata to each dataset?
            xarr.attrs.create('xunits', self.Spectrum.xarr.xunits)
            xarr.attrs.create('xtype', self.Spectrum.xarr.xtype)
            
            if extras is not None:
                for extra in extras:
                    f.create_dataset(extra[0], data = extra[1])
                
            f.close()
            