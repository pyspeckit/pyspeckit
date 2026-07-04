from __future__ import print_function
import re
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


def strip_multispec_wcs(header):
    """
    Remove IRAF echelle/multispec WCS keywords from a header (in-place).

    Spectra loaded from IRAF multispec (echelle) files - e.g., via
    `pyspeckit.wrappers.load_IRAF_multispec` - inherit the original file's
    multispec WCS keywords (``WATi_jjj``, ``CTYPEn = 'MULTISPE'``, etc.).
    When a single order is written back out as a standard linear-WCS
    spectrum, those keywords are stale and cause the file to be
    mis-identified as an echelle spectrum when read back in.
    """
    if 'multispec' not in str(header.get('WAT0_001', '')):
        return header

    for key in list(header.keys()):
        # WATi_jjj: the multispec dispersion specification itself
        # CDi_j / LTMi_j: placeholder (identity) transforms in multispec
        # files; if left in place they override the linear CRVAL/CDELT
        # keywords when the file is read back in
        if key and re.match(r'^(WAT[0-9]+_[0-9]+|CD[0-9]+_[0-9]+|LTM[0-9]+_[0-9]+|LTV[0-9]+)$',
                            key):
            del header[key]

    for key in ('WCSDIM', 'WAXMAP01'):
        if key in header:
            del header[key]

    # axis-2+ multispec keywords describe the order stacking, which no
    # longer exists in a single written spectrum
    for ii in range(2, 8):
        if header.get('CTYPE%i' % ii) == 'MULTISPE':
            for prefix in ('CTYPE', 'CUNIT', 'CRVAL', 'CRPIX', 'CDELT'):
                if '%s%i' % (prefix, ii) in header:
                    del header['%s%i' % (prefix, ii)]

    return header


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

        # copy so that writing does not modify the spectrum's own header
        # (multiple spectra - e.g. echelle orders - may share one header)
        header = self.Spectrum.header.copy()
        # remove stale IRAF multispec WCS keywords (e.g. for a single
        # echelle order loaded with pyspeckit.wrappers.load_IRAF_multispec)
        strip_multispec_wcs(header)
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
