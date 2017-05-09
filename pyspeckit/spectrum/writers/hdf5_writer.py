from __future__ import print_function
import os
from . import Writer

class write_hdf5(Writer):

    def write_data(self, filename=None, newsuffix='out', extras=None,
                   overwrite=True):
        """
        Write information to hdf5 file.

        extras = [(dataset name, data), (dataset name, data)]

        To do: leave option of writing to groups (for model components especially?)
        """


        try:
            import h5py
        except ImportError:
            raise ImportError("Cannot write to hdf5 - h5py import failed.")

        if filename is None:
            fn = "{0}_{1}.hdf5".format(self.Spectrum.fileprefix, newsuffix)
        else:
            fn = filename

        if overwrite:
            f = h5py.File(fn, 'w')
        else:
            i = 1
            while os.path.exists(fn):
                fn = "{0}_fit({1}).hdf5".format(self.Spectrum.fileprefix, i)
                i += 1

        # By default, write out x-array, spectrum, errors
        xarr = f.create_dataset('xarr', data=self.Spectrum.xarr)
        data = f.create_dataset('data', data=self.Spectrum.data)
        #error = f.create_dataset('error', data=self.Spectrum.error)

        # Add metadata to each dataset?
        xarr.attrs.create('unit', str(self.Spectrum.xarr.unit.to_string()))
        data.attrs.create('unit', self.Spectrum.unit.to_string().encode('utf-8')
                          if hasattr(self.Spectrum.unit, 'to_string')
                          else self.Spectrum.unit)
        data.attrs.create('type', self.Spectrum.ytype)

        if extras is not None:
            for extra in extras:
                f.create_dataset(extra[0], data=extra[1]) # add try-except

        f.close()
