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
            data = f.create_group('data')
            data.create_dataset('xarr', data = self.Spectrum.xarr)
            data.create_dataset('data', data = self.Spectrum.data)
            data.create_dataset('err', data = self.Spectrum.error)
            
            #if extras is not None:
            #    xtra = f.create_group('extras')
            #    
            #    for extra in extras:
            #        grp.create_dataset('xarr', data = self.Spectrum.xarr[self.Spectrum.specfit.gx1:self.Spectrum.specfit.gx2])
            #    
            #    grp.create_dataset('model_spectrum', data = self.Spectrum.specfit.model)
            #    grp.create_dataset('model_params', data = self.Spectrum.specfit.modelpars)
            #    f.create_dataset('model_errs', data = self.Spectrum.specfit.modelerrs)
            #    f.create_dataset('residuals', data = self.Spectrum.specfit.residuals)
            #    for i, element in enumerate(self.Spectrum.specfit.modelcomponents):
            #        f.create_dataset('model_component{0}'.format(i), data = self.Spectrum.specfit.modelcomponents[i])
                
            f.close()