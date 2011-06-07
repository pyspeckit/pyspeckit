import os
import numpy as np
from . import Writer

fitscheck = True
try: import pyfits
except ImportError: fitscheck = False

class write_fits(Writer):
    def write_data(self, filename = None, newsuffix = 'out', clobber = True, tolerance=1e-8,
            write_error=True ):
        """
        Write spectrum to fits file.
        """
        
        if not fitscheck: 
            print "Cannot write to fits - pyfits import failed."
            return
        
        if filename is None: fn = "{0}_{1}.fits".format(self.Spectrum.fileprefix, newsuffix)
        else: fn = filename

        header = self.Spectrum.header

        if header.get('CD1_1'): del header['CD1_1']

        # Generate a WCS header from the X-array
        if self.Spectrum.xarr._make_header(tolerance=tolerance):
            for k,v in self.Spectrum.xarr.wcshead.iteritems():
                header.update(k,v)
            if write_error:
                data = np.array( [self.Spectrum.data, self.Spectrum.error] )
            else:
                data = self.Spectrum.data
            print "Writing a FITS-standard (linear-x-axis) spectrum to %s" % (fn)
        else:
            # if no header, overwrite header parameters that would be deceptive
            for k,v in self.Spectrum.xarr.wcshead.iteritems():
                if v is None:
                    if header.get(k): del header[k]
                else:
                    header.update(k,v)
            if write_error:
                data = np.array( [self.Spectrum.xarr, self.Spectrum.data, self.Spectrum.error] )
            else:
                data = np.array( [self.Spectrum.xarr, self.Spectrum.data] )
            print "Writing a nonlinear X-axis spectrum to %s (header keywords are not FITS-compliant)" % (fn)

        HDU = pyfits.PrimaryHDU(data=data, header=header)
        
        HDU.verify('fix')
        HDU.writeto(fn,clobber=clobber)
