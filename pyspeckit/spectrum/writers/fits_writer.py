from __future__ import print_function
from six import iteritems
import astropy
import numpy as np
from . import Writer
import pyspeckit
import warnings

fitscheck = True
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
except ImportError:
    fitscheck = False

class write_fits(Writer):
    def write_data(self, filename=None, newsuffix='out', overwrite=True,
                   tolerance=1e-8, write_error=True, **kwargs):
        """
        Write spectrum to fits file.
        """
        
        if not fitscheck:
            raise ImportError("Could not import FITS handler (astropy.io.fits or pyfits).")
        
        if filename is None:
            fn = "{0}_{1}.fits".format(self.Spectrum.fileprefix, newsuffix)
        else:
            fn = filename

        header = self.Spectrum.header
        header['ORIGIN'] = 'pyspeckit version %s' % pyspeckit.__version__
        header['OBJECT'] = self.Spectrum.specname
        unit = self.Spectrum.unit or self.Spectrum.header.get('BUNIT')
        if unit is not None:
            try:
                header['BUNIT'] = unit.to_string()
            except AttributeError:
                header['BUNIT'] = unit

        #header_nowcs = wcs_utils.strip_wcs_from_header(header)
        #header.insert(2, pyfits.Card(keyword='NAXIS', value=1))
        #header.insert(3, pyfits.Card(keyword='NAXIS1', value=len(self.Spectrum)))

        # Generate a WCS header from the X-array
        if self.Spectrum.xarr._make_header(tolerance=tolerance):
            for k,v in iteritems(self.Spectrum.xarr.wcshead):
                if v is not None:
                    try:
                        header[k] = v
                    except ValueError:
                        try:
                            #v is a Quantity
                            header[k] = v.value
                        except AttributeError:
                            #v is a Unit
                            header[k] = v.to_string()

            if write_error:
                data = np.array( [self.Spectrum.data, self.Spectrum.error] )
            else:
                data = self.Spectrum.data
            print(("Writing a FITS-standard (linear-x-axis) spectrum to %s" % (fn)))
        else:
            # if no header, overwrite header parameters that would be deceptive
            for k,v in iteritems(self.Spectrum.xarr.wcshead):
                if v is None:
                    if header.get(k): del header[k]
                else:
                    header[k] = v
            if write_error:
                data = np.array( [self.Spectrum.xarr, self.Spectrum.data, self.Spectrum.error] )
            else:
                data = np.array( [self.Spectrum.xarr, self.Spectrum.data] )
            warnings.warn("Writing a nonlinear X-axis spectrum to %s (header keywords are not FITS-compliant)" % (fn))

        try:
            HDU = pyfits.PrimaryHDU(data=data, header=header)
        except AttributeError:
            print("Strange header error. Attempting workaround.")
            HDU = pyfits.PrimaryHDU(data=data,
                                    header=pyfits.Header([pyfits.card.Card(k,v)
                                                          for k,v in
                                                          iteritems(header)]))
        
        HDU.verify('fix')
        if astropy.version.major >= 2 or (astropy.version.major==1 and astropy.version.minor>=3):
            HDU.writeto(fn, overwrite=overwrite, output_verify='fix', **kwargs)
        else:
            HDU.writeto(fn, clobber=overwrite, output_verify='fix', **kwargs)
