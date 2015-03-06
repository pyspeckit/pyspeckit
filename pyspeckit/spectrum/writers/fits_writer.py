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
    def write_data(self, filename=None, newsuffix='out', clobber=True,
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
            header['BUNIT'] = unit

        #header_nowcs = wcs_utils.strip_wcs_from_header(header)
        #header.insert(2, pyfits.Card(keyword='NAXIS', value=1))
        #header.insert(3, pyfits.Card(keyword='NAXIS1', value=len(self.Spectrum)))

        # Generate a WCS header from the X-array
        if self.Spectrum.xarr._make_header(tolerance=tolerance):
            # print('iteratig over:', self.Spectrum.xarr.wcshead)
            for k,v in self.Spectrum.xarr.wcshead.iteritems():
                # print('v is type: %s and has a value: %s' % (type(v), v))
                if v is not None:
                    try:
                        header[k] = v
                    except ValueError:
                        try:
                            header[k] = v.value
                        except AttributeError:
                            pass
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
                    header[k] = v
            if write_error:
                data = np.array( [self.Spectrum.xarr, self.Spectrum.data, self.Spectrum.error] )
            else:
                data = np.array( [self.Spectrum.xarr, self.Spectrum.data] )
            warnings.warn("Writing a nonlinear X-axis spectrum to %s (header keywords are not FITS-compliant)" % (fn))

        try:
            HDU = pyfits.PrimaryHDU(data=data, header=header)
        except AttributeError:
            print "Strange header error. Attempting workaround."
            HDU = pyfits.PrimaryHDU(data=data, header=pyfits.Header([pyfits.card.Card(k,v)  for k,v in header.iteritems()]))
        
        HDU.verify('fix')
        HDU.writeto(fn, clobber=clobber, output_verify='fix', **kwargs)
