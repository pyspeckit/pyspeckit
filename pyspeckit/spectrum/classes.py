"""
==========================
PySpecKit Spectrum classes
==========================

The spectrum module consists of the Spectrum class, with child classes ObsBlock
and Spectra for multi-spectrum analysis of different types.

The Spectrum class is the main functional code.
ObsBlocks are containers of multiple spectra of different objects
The Spectra class is a container of multiple spectra of the *same* object at
different wavelengths/frequencies

.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>
.. moduleauthor:: Jordan Mirocha <mirochaj@gmail.com>
"""
from __future__ import print_function
import numpy as np
from six import iteritems
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
from . import smooth as sm
from . import readers
from . import plotters
from . import writers
from . import baseline
from . import units
from . import measurements
from . import speclines
from . import interpolation
from . import moments as moments_module
from . import fitters
from . import history
import copy
from astropy import log
from ..specwarnings import warn, PyspeckitWarning
try:
    import atpy
    atpyOK = True
except ImportError:
    atpyOK = False

# specutils -> legacy specutils
# try:
#     import specutils
#     specutilsOK = True
# except ImportError:
#     specutilsOK = False


try:
    from Spectrum1D import Spectrum1D # inherit from astropy
except ImportError:
    Spectrum1D = object

try:
    import astropy.units as u
except ImportError:
    u = None

class BaseSpectrum(object):

    from .interpolation import interpnans

    def __init__(self, filename_or_magic=None, filetype=None, xarr=None,
                 data=None, error=None, header=None, doplot=False,
                 maskdata=True, unit=None, plotkwargs={}, xarrkwargs={},
                 model_registry=None, filename=None,
                 **kwargs):
        """
        Create a Spectrum object.

        Must either pass in a filename or ALL of xarr, data, and header, plus
        optionally error.

        kwargs are passed to the file reader

        Parameters
        ----------
        filename_or_magic : string or something else
            The filename or something with an hdu attribute.  If data, xarr, and error are
            specified, leave filename blank.
        filetype : string
            Specify the file type (only needed if it cannot be automatically
            determined from the filename)
        xarr : `units.SpectroscopicAxis` or `np.ndarray`
            The X-axis of the data.  If it is an np.ndarray, you must pass
            `xarrkwargs` or a valid header if you want to use any of the unit
            functionality.
        data : `np.ndarray`
            The data array (must have same length as xarr)
        error : `np.ndarray`
            The error array (must have same length as the data and xarr arrays)
        header : `pyfits.Header` or dict
            The header from which to read unit information.  Needs to be a
            `pyfits.Header` instance or another dictionary-like object with the
            appropriate information
        maskdata : boolean
            turn the array into a masked array with all nan and inf values masked
        doplot : boolean
            Plot the spectrum after loading it?
        plotkwargs : dict
            keyword arguments to pass to the plotter
        xarrkwargs : dict
            keyword arguments to pass to the SpectroscopicAxis initialization
            (can be used in place of a header)
        unit : str
            The data unit
        filename : string
            The file to read the spectrum from.  If data, xarr, and error are
            specified, leave filename blank.

        Examples
        --------

        >>> sp = pyspeckit.Spectrum(data=np.random.randn(100),
                    xarr=np.linspace(-50, 50, 100), error=np.ones(100)*0.1,
                    xarrkwargs={'unit':'km/s', 'refX':4.829, 'refX_unit':'GHz',
                        'xtype':'VLSR-RAD'}, header={})

        >>> xarr = pyspeckit.units.SpectroscopicAxis(np.linspace(-50,50,100),
                    units='km/s', refX=6562.83, refX_unit='angstroms')
        >>> data = np.random.randn(100)*5 + np.random.rand(100)*100
        >>> err = np.sqrt(data/5.)*5. # Poisson noise
        >>> sp = pyspeckit.Spectrum(data=data, error=err, xarr=xarr, header={})

        >>> # if you already have a simple fits file
        >>> sp = pyspeckit.Spectrum('test.fits')
        """
        if filename_or_magic is not None:
            if hasattr(filename_or_magic, 'hdu'):
                return self.from_hdu(filename_or_magic.hdu)
            elif filename is None:
                filename = filename_or_magic
            else:
                raise ValueError("filename_or_magic was specified incorrectly")
        if filename:
            if error is not None:
                raise ValueError("When reading from a file, you cannot specify"
                                 "the error as an array.  Instead, set it "
                                 "separately after reading the file, e.g.: \n"
                                 "sp = Spectrum(filename)\n"
                                 "sp.error[:] = rms")
            if xarr is not None:
                raise ValueError("Cannot specify xarr when reading from a "
                                 "file.  If the xarr in the file is incorrect,"
                                 "change it after reading the file in, i.e., "
                                 "set sp.xarr on another line.")

            if filetype is None:
                suffix = filename.rsplit('.',1)[1]
                if suffix in readers.suffix_types:
                    # use the default reader for that suffix
                    filetype = readers.suffix_types[suffix][0]
                    reader = readers.readers[filetype]
                else:
                    raise TypeError("File with suffix %s is not recognized." % suffix)
            else:
                if filetype in readers.readers:
                    reader = readers.readers[filetype]
                else:
                    raise TypeError("Filetype %s not recognized" % filetype)

            self.data,self.error,self.xarr,self.header = reader(filename,**kwargs)

            # these should probably be replaced with registerable function s...
            if filetype in ('fits','tspec','pyfits','sdss'):
                self.parse_header(self.header)
            elif filetype == 'txt':
                self.parse_text_header(self.header)
            elif filetype in ('hdf5', 'h5'):
                self.parse_hdf5_header(self.header)

            if isinstance(filename,str):
                # Everything prior to .fits or .txt
                self.fileprefix = filename.rsplit('.', 1)[0]
        elif xarr is not None and data is not None:
            # technically, this is unpythonic.  But I don't want to search for
            # all 10 attributes required.
            if issubclass(type(xarr),units.SpectroscopicAxis):
                self.xarr = xarr
            else:
                self.xarr = units.SpectroscopicAxis(xarr, **xarrkwargs)
            self.data = data
            if error is not None:
                self.error = error
            else:
                self.error = np.zeros_like(data)
            if hasattr(header,'get'):
                if not isinstance(header, pyfits.Header):
                    cards = [pyfits.Card(k, header[k]) for k in header]
                    self.header = pyfits.Header(cards)
                else:
                    self.header = header
            else: # set as blank
                warn("WARNING: No header given.  Creating an empty one.",
                     PyspeckitWarning)
                self.header = pyfits.Header()
            self.parse_header(self.header)
        else:
            raise ValueError("Must either give a filename or xarr and data "
                             "keywords to instantiate a pyspeckit.Spectrum")

        if hasattr(self.data,'unit'):
            # TODO: use the quantity more appropriately
            self.unit = str(self.data.unit)
            self.data = self.data.value

        if hasattr(self.error, 'unit'):
            # errors have to have the same units as data, but then should be
            # converted to arrays.  Long term, we'd like to have everything
            # be treated internally as a Quantity, but... not yet.
            self.error = self.error.to(self.unit).value

        if maskdata:
            if hasattr(self.data,'mask'):
                self.data.mask += np.isnan(self.data) + np.isinf(self.data)
                if hasattr(self.error,'mask'):
                    self.error.mask += np.isnan(self.data) + np.isinf(self.data)
            else:
                self.data = np.ma.masked_where(np.isnan(self.data) + np.isinf(self.data), self.data)
                self.error = np.ma.masked_where(np.isnan(self.data) + np.isinf(self.data), self.error)

        # it is very important that this be done BEFORE the spectofit is set!
        self._sort()
        self.plotter = plotters.Plotter(self)
        self._register_fitters(registry=model_registry)
        self.specfit = fitters.Specfit(self, Registry=self.Registry)
        self.baseline = baseline.Baseline(self)
        self.speclines = speclines

        # Special.  This needs to be modified to be more flexible; for now I need it to work for nh3
        self.plot_special = None
        self.plot_special_kwargs = {}

        if unit is not None:
            self._unit = unit
        elif not hasattr(self, '_unit'):
            self._unit = u.dimensionless_unscaled

        if doplot:
            self.plotter(**plotkwargs)

    @property
    def unit(self):
        return self._unit

    @property
    def units(self):
        log.warning("'units' is deprecated; please use 'unit'", DeprecationWarning)
        return self._unit

    @unit.setter
    def unit(self, value):
        self._unit = value

    @units.setter
    def units(self, value):
        log.warning("'units' is deprecated; please use 'unit'", DeprecationWarning)
        self._unit = value

    @property
    def data_quantity(self):
        return u.Quantity(self.data, unit=self.unit)

    @data_quantity.setter
    def data_quantity(self, value):
        if not hasattr(value, 'unit'):
            raise ValueError("To set the data to a Quantity value, it must "
                             "have a unit.")
        if hasattr(self.data, 'mask') and not hasattr(value, 'mask'):
            raise ValueError("The original data had a mask.  You must use "
                             "a masked array to set the data value.")
        self.data = value.value
        self.unit = value.unit

    @property
    def error_quantity(self):
        return u.Quantity(self.error, unit=self.unit)

    def _register_fitters(self, registry=None):
        """
        Register fitters independently for each spectrum instance

        This approach allows you to add fitters to a given Spectrum instance
        without modifying the default registry
        """
        self.Registry = fitters.Registry()
        if registry is None:
            registry = fitters.default_Registry
        elif not isinstance(registry, fitters.Registry):
            raise TypeError("registry must be an instance of the fitters.Registry class")

        for modelname, model in iteritems(registry.multifitters):
            self.Registry.add_fitter(modelname, model,
                                     registry.npars[modelname],
                                     key=registry.associated_keys.get(modelname))

    def _sort(self):
        """
        Make sure X axis is monotonic.
        """
        if self.xarr.dxarr.min() < 0:
            argsort = np.argsort(self.xarr)
            self.data = self.data[argsort]
            self.error = self.error[argsort]
            self.xarr = self.xarr[argsort]
            self.xarr.make_dxarr()

    def write(self,filename,type=None,**kwargs):
        """
        Write the spectrum to a file.  The available file types are listed
        in spectrum.writers.writers

        type - what type of file to write to?  If not specified, will attempt
        to determine type from suffix
        """
        if type:
            self.writer = writers.writers[type](self)
        else:
            suffix = filename.rsplit('.',1)[1]
            if suffix in writers.suffix_types:
                # use the default reader for that suffix
                filetype = writers.suffix_types[suffix][0]
                self.writer = writers.writers[filetype](self)
            else:
                raise TypeError("File with suffix %s is not recognized." % suffix)
        self.writer(filename=filename,**kwargs)

    def parse_text_header(self,Table):
        """
        Grab relevant parameters from a table header (xaxis type, etc)

        This function should only exist for Spectrum objects created from
        .txt or other atpy table type objects
        """
        self.Table = Table

        xtype = Table.data.dtype.names[Table.xaxcol]
        if xtype in units.xtype_dict.values():
            self.xarr.xtype = xtype
            unit = Table.columns[xtype].unit
            self.xarr.set_unit(unit)
        elif xtype in units.xtype_dict:
            self.xarr.xtype = units.xtype_dict[xtype]
            unit = Table.columns[xtype].unit
            self.xarr.set_unit(unit)
        else:
            warn("Warning: Invalid xtype in text header - this may mean no "
                 "text header was available.  X-axis units will be pixels "
                 "unless you set them manually "
                 "(e.g., sp.xarr=SpectroscopicAxis(sp.xarr.value, unit='angstroms')")
            self.xarr.xtype = 'pixels'
            self.xarr.set_unit(u.pixel)
            #raise ValueError("Invalid xtype in text header")
        self.ytype = Table.data.dtype.names[Table.datacol]
        try:
            self.unit = Table.columns[self.ytype].unit
        except ValueError:
            self.unit = None
            pass # astropy 0.2.dev11 introduces an incompatibility here
        self.header = pyfits.Header()
        self._update_header()

    def _update_header(self):
        self.header['CUNIT1'] = self.xarr.unit.to_string()
        self.header['CTYPE1'] = self.xarr.xtype
        self.header['BUNIT'] = self.unit
        self.header['BTYPE'] = self.ytype

    def parse_hdf5_header(self, hdr):
        """
        HDF5 reader will create a hdr dictionary from HDF5 dataset attributes
        if they exist.  This routine will convert that dict to a pyfits header
        instance.

        .. todo:: Move this to the hdf5 reader?
        """

        self.xarr.xtype = hdr['xtype']
        self.xarr.xunit = hdr['xunit']
        self.ytype = hdr['ytype']
        self.unit = hdr['yunit']
        self.header = pyfits.Header()
        self.header['CUNIT1'] = self.xarr.xunit
        self.header['CTYPE1'] = self.xarr.xtype
        self.header['BUNIT'] = self.ytype
        self.header['BTYPE'] = self.unit

    def parse_header(self,hdr,specname=None):
        """
        Parse parameters from a .fits header into required spectrum structure
        parameters

        .. todo:: This should be moved to the FITSSpectrum subclass when that is available
        """

        if hdr.get('BUNIT'):
            try:
                self.unit = u.Unit(hdr.get('BUNIT').strip())
            except (u.UnitsError, ValueError):
                self.unit = hdr.get('BUNIT').strip()
        elif not hasattr(self, 'unit') or (hasattr(self,'unit') and self.unit
                                           is None):
            self.unit = 'undefined'

        if hdr.get('BTYPE'):
            self.ytype = hdr.get('BTYPE').strip()
        else:
            self.ytype = 'data'

        if specname is not None:
            self.specname = specname
        elif hdr.get('OBJECT'):
            self.specname = hdr.get('OBJECT')
        elif hdr.get('OBJNAME'):
            self.specname = hdr.get('OBJNAME')
        else:
            self.specname = ''

    def measure(self, z=None, d=None, fluxnorm=None, miscline=None,
                misctol=10.0, ignore=None, derive=True, **kwargs):
        """
        Initialize the measurements class - only do this after you have run a
        fitter otherwise pyspeckit will be angry!
        """
        self.measurements=measurements.Measurements(self, z=z, d=d,
                                                    fluxnorm=fluxnorm,
                                                    miscline=miscline,
                                                    misctol=misctol,
                                                    ignore=ignore,
                                                    derive=derive, **kwargs)

    def crop(self, x1, x2, unit=None, **kwargs):
        """
        Replace the current spectrum with a subset from x1 to x2 in current
        units

        Fixes CRPIX1 and baseline and model spectra to match cropped data spectrum

        """

        # do slice (this code is redundant... need to figure out how to fix that)
        if not hasattr(x1, 'unit') and unit is not None:
            x1pix = np.argmin(np.abs(x1-self.xarr.as_unit(unit).value))
            x2pix = np.argmin(np.abs(x2-self.xarr.as_unit(unit).value))
        elif hasattr(x1, 'unit') and unit is not None:
            raise ValueError("If you give x1,x2 as quantities, don't specify "
                             "the X-axis unit (it must be equivalent, though).")
        elif not hasattr(x1, 'unit') and unit is None:
            x1pix = x1
            x2pix = x2
        else:
            # Hack: something about the inheritance of xarr prevents equivalent
            # unit arithmetic
            x1pix = np.argmin(np.abs(x1-u.Quantity(self.xarr)))
            x2pix = np.argmin(np.abs(x2-u.Quantity(self.xarr)))
            x1 = x1.value
            x2 = x2.value

        if x1pix > x2pix:
            x1pix,x2pix = x2pix,x1pix
        elif x1pix == x2pix:
            raise IndexError("ERROR: Trying to crop to zero size.")

        self = self.slice(x1pix, x2pix, unit='pixels', copy=False, xcopy=True,
                          **kwargs)
        # a baseline spectrum is always defined, even if it is all zeros
        # this is needed to prevent size mismatches.  There may be a more
        # elegant way to do this...
        # this is taken care of in Slice now self.baseline.crop(x1pix,x2pix)
        # this is taken care of in Slice now self.specfit.crop(x1pix,x2pix)
        if hasattr(self.specfit, 'fitter'):
            self.specfit._full_model()

        if hasattr(self,'header'):
            history.write_history(self.header,
                                  "CROP: Cropped from %g to %g (pixel %i to %i)"
                                  % (x1,x2,x1pix,x2pix))

            if self.header.get('CRPIX1'):
                self.header['CRPIX1'] = self.header.get('CRPIX1') - x1pix
                history.write_history(self.header,
                                      "CROP: Changed CRPIX1 from %f to %f"
                                      % (self.header.get('CRPIX1')+x1pix,self.header.get('CRPIX1')))

    def slice(self, start=None, stop=None, unit='pixel', copy=True, xcopy=True,
              preserve_fits=False):
        """
        Slicing the spectrum

        .. WARNING:: this is the same as cropping right now, but it returns a
            copy instead of cropping inplace

        Parameters
        ----------
        start : numpy.float or int or astropy quantity
            start of slice
        stop :  numpy.float or int or astropy quantity
            stop of slice
        unit : str
            allowed values are any supported physical unit, 'pixel'
        copy : bool
            Return a 'view' of the data or a copy?
        preserve_fits : bool
            Save the fitted parameters from self.fitter?
        """

        if hasattr(start, 'unit'):
            start_ind = self.xarr.x_to_pix(start)
        elif unit in ('pixel','pixels'):
            start_ind = start
        else:
            start_ind = self.xarr.x_to_pix(start, xval_unit=unit)

        if hasattr(stop, 'unit'):
            stop_ind = self.xarr.x_to_pix(stop)
        elif unit in ('pixel','pixels'):
            stop_ind = stop
        else:
            stop_ind = self.xarr.x_to_pix(stop, xval_unit=unit)

        if start_ind > stop_ind:
            start_ind,stop_ind = stop_ind,start_ind

        spectrum_slice = slice(start_ind,stop_ind)

        log.debug("Slicing from {start} to {stop} with unit {unit} and copy="
                  "{copy}, xcopy={xcopy}, preserve_fits={preserve_fits}."
                  "start_ind = {start_ind}, stop_ind= {stop_ind}"
                  .format(start=start, stop=stop, unit=unit, copy=copy,
                          xcopy=xcopy, preserve_fits=preserve_fits,
                          start_ind=start_ind, stop_ind=stop_ind))

        if copy:
            sp = self.copy()
        else:
            sp = self
        sp.data = sp.data[spectrum_slice]
        if sp.error is not None:
            sp.error = sp.error[spectrum_slice]
        if copy or xcopy:
            sp.xarr = sp.xarr[spectrum_slice].copy()
        else:
            sp.xarr = sp.xarr[spectrum_slice]

        if copy:
            # create new specfit / baseline instances (otherwise they'll be the wrong length)
            sp._register_fitters(registry=self.Registry)
            sp.baseline = baseline.Baseline(sp)
            sp.specfit = fitters.Specfit(sp,Registry=sp.Registry)
        else:
            # inplace modification
            sp.baseline.crop(start_ind, stop_ind)
            sp.specfit.crop(start_ind, stop_ind)

        if preserve_fits:
            sp.specfit.modelpars = self.specfit.modelpars
            sp.specfit.parinfo = self.specfit.parinfo
            sp.baseline.baselinepars = self.baseline.baselinepars
            sp.baseline.order = self.baseline.order


        return sp

    # For Spectrum1D compatibility, flux = data
    @property
    def flux(self):
        """
        The data in the spectrum (flux = data, for compatibility with astropy's
        Spectrum1D object).
        """
        return self.data

    @flux.setter
    def flux(self,value):
        self.data = value

    def __getitem__(self, indx):
        """
        Slice the data using pixel units (not quite the same as self.slice
        because indx can be a slice object or numbers)
        """

        sp = copy.copy(self)
        sp.data = sp.data.__getitem__(indx)
        if sp.error is not None:
            sp.error = sp.error.__getitem__(indx)
        sp.xarr = sp.xarr.__getitem__(indx)

        # this should be done by deepcopy, but deepcopy fails with current pyfits
        sp.plotter = self.plotter.copy(parent=sp)
        sp.plotter.Spectrum = sp
        sp.specfit = self.specfit.copy(parent=sp, registry=sp.Registry)
        sp.specfit.Spectrum = sp
        sp.specfit.Spectrum.plotter = sp.plotter
        sp.baseline = self.baseline.copy(parent=sp)
        sp.baseline.Spectrum = sp
        sp.baseline.Spectrum.plotter = sp.plotter

        return sp

    def downsample(self, dsfactor):
        """
        Downsample the spectrum (and all of its subsidiaries) without smoothing

        Parameters
        ----------
        dsfactor : int
            Downsampling Factor
        """

        dsfactor = round(dsfactor)  # ensure int

        self.data = self.data[::dsfactor]
        self.xarr = self.xarr[::dsfactor]
        if len(self.xarr) != len(self.data):
            raise ValueError("Convolution resulted in different X and Y array lengths.  Convmode should be 'same'.")
        if self.error is not None:
            self.error = self.error[::dsfactor]
        self.baseline.downsample(dsfactor)
        self.specfit.downsample(dsfactor)

    def smooth(self,smooth,downsample=True,**kwargs):
        """
        Smooth the spectrum by factor `smooth`.


        Documentation from the :mod:`smooth` module:

        """
        smooth = int(round(smooth))
        self.data = sm.smooth(self.data,smooth,downsample=downsample,**kwargs)

        if downsample:
            self.xarr = self.xarr[::smooth]
            if len(self.xarr) != len(self.data):
                raise ValueError("Convolution resulted in different X and Y array lengths.  Convmode should be 'same'.")
            if self.error is not None:
                self.error = sm.smooth(self.error,smooth,**kwargs)
            self.baseline.downsample(smooth)
            self.specfit.downsample(smooth)

            self._smooth_header(smooth)
    smooth.__doc__ += "sm.smooth doc: \n" + sm.smooth.__doc__

    def _smooth_header(self,smooth):
        """
        Internal - correct the FITS header parameters when smoothing
        """
        if self.header.get('CDELT1') is not None and self.header.get('CRPIX1') is not None:
            self.header['CDELT1'] = self.header.get('CDELT1') * float(smooth)
            self.header['CRPIX1'] = self.header.get('CRPIX1') / float(smooth)

            history.write_history(self.header,"SMOOTH: Smoothed and downsampled spectrum by factor %i" % (smooth))
            history.write_history(self.header,"SMOOTH: Changed CRPIX1 from %f to %f" % (self.header.get('CRPIX1')*float(smooth),self.header.get('CRPIX1')))
            history.write_history(self.header,"SMOOTH: Changed CDELT1 from %f to %f" % (self.header.get('CRPIX1')/float(smooth),self.header.get('CRPIX1')))

    def _shape(self):
        """
        Return the data shape (a property of the Spectrum)
        """
        return self.data.shape

    shape = property(_shape)

    def __len__(self):
        return len(self.data)

    def __repr__(self):
        if hasattr(self,'specname'):
            name = " named %s" % self.specname
        else:
            name = ""
        return (r'<%s object%s over spectral range %6.5g : %6.5g %s and flux range = [%2.1f, %2.1f] %s at %s>' %
                (self.__class_name__, name, self.xarr.min().value, self.xarr.max().value,
                 self.xarr.unit, self.data.min(), self.data.max(), self.unit,
                 str(hex(self.__hash__()))))


    def copy(self,deep=True):
        """
        Create a copy of the spectrum with its own plotter, fitter, etc.
        Useful for, e.g., comparing smoothed to unsmoothed data
        """

        newspec = copy.copy(self)
        if deep:
            newspec.xarr = copy.copy(self.xarr)
            newspec.data = copy.copy(self.data)
            if self.error is not None:
                newspec.error = copy.copy(self.error)

        newspec.header = copy.copy(self.header)
        newspec.plotter = self.plotter.copy(parent=newspec)
        newspec._register_fitters(registry=self.Registry)
        newspec.specfit = self.specfit.copy(parent=newspec, registry=newspec.Registry)
        newspec.specfit.Spectrum.plotter = newspec.plotter
        newspec.baseline = self.baseline.copy(parent=newspec)
        newspec.baseline.Spectrum.plotter = newspec.plotter

        return newspec

    def stats(self, statrange=(), interactive=False):
        """
        Return some statistical measures in a dictionary (somewhat self-explanatory)

        Parameters
        ----------
        statrange : 2-element tuple
            X-range over which to perform measures
        interactive : bool
            specify range interactively in plotter
        """

        if len(statrange) == 2:
            pix1 = self.xarr.x_to_pix(statrange[0])
            pix2 = self.xarr.x_to_pix(statrange[1])
            if pix1 > pix2:
                pix1,pix2 = pix2,pix1
            elif pix1 == pix2:
                raise ValueError("Invalid statistics range - includes 0 pixels")
            data = self.data[pix1:pix2]
        elif interactive:
            raise NotImplementedError('Not implemented yet.  Probably need to move the stats command into a different module first')
        else:
            data = self.data

        stats = {
            "npts": data.shape[0],
            "std": data.std(),
            "mean": data.mean(),
            "median": np.median(data),
            "min": data.min(),
            "max": data.max(),}
        return stats

    def getlines(self, linetype='radio', **kwargs):
        """
        Access a registered database of spectral lines.  Will add an attribute
        with the name linetype, which then has properties defined by the
        speclines module (most likely, a table and a "show" function to display
        the lines)
        """

        # this is somewhat unreadable, but in short:
        # self.__dict__ is a dictionary that defines the class attributes
        # so, for example, if linetype is radio, this reads:
        # self.radio = speclines.radio.radio_lines(self)
        # or optical:
        # self.optical = speclines.optical.optical_lines(self)
        if linetype not in self.__dict__: # don't replace it if it already exists
            self.__dict__[linetype] = speclines.__dict__[linetype].__dict__[linetype+"_lines"](self,**kwargs)

    def moments(self, unit='km/s', **kwargs):
        """
        Return the moments of the spectrum.  In order to assure that the 1st
        and 2nd moments are meaningful, a 'default' unit is set.  If unit is
        not set, will use current unit.

        *Documentation imported from the moments module:*
        """

        if unit is False or unit is None:
            return moments_module.moments(self.xarr, self.data, **kwargs)
        else:
            return moments_module.moments(self.xarr.as_unit(unit), self.data, **kwargs)

    moments.__doc__ += moments_module.moments.__doc__

    def _operation_wrapper(operation):
        """
        Perform an operation (addition, subtraction, mutiplication, division, etc.)
        after checking for shape matching
        """


        def ofunc(self, other):
            if np.isscalar(other):
                newspec = self.copy()
                newspec.data = operation(newspec.data, other)
                return newspec

            elif hasattr(self,'xarr') and hasattr(other,'xarr'): # purely for readability

                if self._arithmetic_threshold == 'exact':
                    xarrcheck = np.all((self.xarr == other.xarr))
                else:
                    if self._arithmetic_threshold_units is None:
                        # not sure this should ever be allowed
                        xarrcheck = np.all(np.abs(self.xarr-other.xarr).value < self._arithmetic_threshold)
                    else:
                        xarr_u = self.xarr.as_unit(self._arithmetic_threshold_units)
                        other_xarr_u = other.xarr.as_unit(self._arithmetic_threshold_units)
                        xarrcheck = np.all(np.abs(xarr_u - other_xarr_u).value <
                                           self._arithmetic_threshold)

                if self.shape == other.shape and xarrcheck:
                    newspec = self.copy()
                    newspec.data = operation(newspec.data, other.data)
                    return newspec
                elif self.shape != other.shape:
                    raise ValueError("Shape mismatch in data")
                elif not xarrcheck:
                    raise ValueError("X-axes do not match.")
                else:
                    raise Exception("Unexpected Error")
            elif hasattr(self,'shape') and hasattr(other,'shape'):
                # allow array subtraction
                if self.shape != other.shape:
                    raise ValueError("Shape mismatch in data")
                elif hasattr(self, 'xarr'):
                    newspec = self.copy()
                    newspec.data = operation(newspec.data, other)
                    return newspec
                elif hasattr(other, 'xarr'): # is this even possible?
                    newspec = other.copy()
                    newspec.data = operation(self, other.data)
                    return newspec
                else:
                    raise ValueError("Data shapes match but somehow neither side is a Spectrum")
            else:
                raise ValueError("Data types are incompatible")

        return ofunc

    @property
    def _arithmetic_threshold(self):
        return self._arithmetic_threshold_value

    @_arithmetic_threshold.setter
    def _arithmetic_threshold(self, value, unit=None):
        self._arithmetic_threshold_value = value
        if unit is None:
            self._arithmetic_threshold_units = self.xarr.unit
        else:
            self._arithmetic_threshold_units = unit

    _arithmetic_threshold_value = 'exact'
    _arithmetic_threshold_units = None

    __add__ = _operation_wrapper(np.add)
    __radd__ = _operation_wrapper(np.add)
    __sub__ = _operation_wrapper(np.subtract)
    __mul__ = _operation_wrapper(np.multiply)
    __div__ = _operation_wrapper(np.divide)

class SingleSpectrum(BaseSpectrum):
    __class_name__ = 'SingleSpectrum'

    @classmethod
    def from_spectrum1d(cls, spec1d):
        """
        Tool to load a pyspeckit Spectrum from a specutils object

        Examples
        --------
        >>> # grab many spectra from a multiextension FITS file
        >>> spectra = specutils.io.fits.read_fits_spectrum1d('AAO.fits')
        >>> sp = pyspeckit.Spectrum.from_spectrum1d(spectra[0])

        >>> # open a single spectrum that could have been opened directly with pyspeckit
        >>> spectrum = specutils.io.fits.read_fits_spectrum1d('gbt_1d.fits')
        >>> sp = pyspeckit.Spectrum.from_spectrum1d(spectrum)
        """
        xarr = units.SpectroscopicAxis(spec1d.spectral_axis)

        data = spec1d.flux
        error = spec1d.uncertainty

        return cls(data=data, error=error, xarr=xarr,
                   unit=spec1d._unit, header=pyfits.Header())


    @classmethod
    def from_hdu(cls, hdu):
        """
        Create a pyspeckit Spectrum object from an HDU
        """

        spec,errspec,XAxis,hdr = readers.open_1d_pyfits(hdu)
        return cls(data=spec, error=errspec, xarr=XAxis, header=hdr)

    __class_name__ = 'BaseSpectrum'


class Spectrum(SingleSpectrum):
    """
    The core class for the spectroscopic toolkit.  Contains the data and error
    arrays along with wavelength / frequency / velocity information in various
    formats.
    """
    # just does what BaseSpectrum and SingleSpectrum do
    __class_name__ = 'Spectrum'


class Spectra(BaseSpectrum):
    """
    A list of individual Spectrum objects.  Intended to be used for
    concatenating different wavelength observations of the SAME OBJECT.  Can be
    operated on just like any Spectrum object, incuding fitting.  Useful for
    fitting multiple lines on non-continguous axes simultaneously.  Be wary of
    plotting these though...

    Can be indexed like python lists.

    X array is forcibly sorted in increasing order
    """
    __class_name__ = 'Spectra'

    def __init__(self, speclist, xunit=None, model_registry=None, **kwargs):
        """
        """
        log.info("Creating spectra")
        speclist = list(speclist)
        for ii,spec in enumerate(speclist):
            if type(spec) is str:
                spec = Spectrum(spec)
                speclist[ii] = spec

        self.speclist = speclist

        if xunit is None:
            xunit = speclist[0].xarr.unit
        else:
            xunit = u.Unit(xunit)

        if xunit is not None and xunit.is_equivalent(u.km/u.s):
            refX = speclist[0].xarr.refX
            if refX is None:
                warn("Combining spectra with velocity coordinates, "
                     "but refX is None")
            for spec in speclist[1:]:
                if spec.xarr.refX != refX:
                    raise ValueError("When combining spectra in velocity coordinates, "
                                     "they must have the same reference frequency.")

        log.info("Concatenating data")

        self.xarr = units.SpectroscopicAxes([sp.xarr.as_unit(xunit) for sp in speclist])
        self.xarr.set_unit(u.Unit(xunit))
        self.xarr.xtype = u.Unit(xunit)
        self.data = np.ma.concatenate([sp.data for sp in speclist])
        self.error = np.ma.concatenate([sp.error for sp in speclist])
        self._sort()

        self.header = pyfits.Header()
        for spec in speclist:
            for key,value in spec.header.items():
                try:
                    self.header[key] = value
                except (ValueError, KeyError):
                    warn("Could not update header KEY=%s to VALUE=%s" % (key,value))

        self.plotter = plotters.Plotter(self)
        self._register_fitters(registry=model_registry)
        self.specfit = fitters.Specfit(self,Registry=self.Registry)
        self.baseline = baseline.Baseline(self)

        self.unit = speclist[0].unit
        for spec in speclist:
            if spec.unit != self.unit:
                raise ValueError("Mismatched unit")

        # Special.  This needs to be modified to be more flexible; for now I need it to work for nh3
        self.plot_special = None
        self.plot_special_kwargs = {}

    @classmethod
    def from_spectrum1d_list(cls, lst):
        """
        Tool to load a collection of pyspeckit Spectra from a specutils list of
        Spectrum1D objects

        Examples
        --------
        >>> # grab many spectra from a multiextension FITS file
        >>> spectra = specutils.io.fits.read_fits_spectrum1d('AAO.fits')
        >>> sp = pyspeckit.Spectrum.from_spectrum1d_list(spectra)
        """
        return cls([Spectrum.from_spectrum1d(spec) for spec in lst])

    def __add__(self,other):
        """
        Adding "Spectra" together will concatenate them
        * WARNING * this will probably fail now that I've implemented direct __add__...
        """
        if type(other) is Spectra or Spectrum:
            self.speclist += other.speclist
        elif type(other) is Spectrum:
            self.speclist += [other.speclist]
        if other.unit != self.unit:
            raise ValueError("Mismatched unit")

        if other.xarr.unit != self.xarr.unit:
            # convert all inputs to same unit
            other.xarr.convert_to_unit(self.xarr.unit)
        self.xarr = units.SpectroscopicAxes([self.xarr,other.xarr])
        self.data = np.concatenate([self.data,other.data])
        self.error = np.concatenate([self.error,other.error])
        self._sort()

    def __getitem__(self,index):
        """
        Can index Spectra to get the component Spectrum objects
        """
        return self.speclist[index]

    def __len__(self):
        """
        len(spectra) != len(spectrum) !
        """
        return len(self.speclist)

    def smooth(self,smooth,**kwargs):
        """
        Smooth the spectrum by factor "smooth".  Options are defined in sm.smooth

        because 'Spectra' does not have a header attribute, don't do anything to it...
        """
        smooth = round(smooth)
        self.data = sm.smooth(self.data,smooth,**kwargs)
        self.xarr = self.xarr[::smooth]
        if len(self.xarr) != len(self.data):
            raise ValueError("Convolution resulted in different X and Y array lengths.  Convmode should be 'same'.")
        self.error = sm.smooth(self.error,smooth,**kwargs)
        self.baseline.downsample(smooth)
        self.specfit.downsample(smooth)

    def fiteach(self,**kwargs):
        """
        Fit each spectrum within the Spectra object
        """
        for sp in self.speclist:
            sp.specfit(**kwargs)

        if atpyOK:
            self.fittable = atpy.Table()
            self.fittable.add_column('name',[sp.specname for sp in self.speclist])
            self.fittable.add_column('amplitude',[sp.specfit.modelpars[0] for sp in self.speclist],unit=self.unit)
            self.fittable.add_column('center',[sp.specfit.modelpars[1] for sp in self.speclist],unit=self.xarr.unit)
            self.fittable.add_column('width',[sp.specfit.modelpars[2] for sp in self.speclist],unit=self.xarr.unit)
            self.fittable.add_column('amplitudeerr',[sp.specfit.modelerrs[0] for sp in self.speclist],unit=self.unit)
            self.fittable.add_column('centererr',[sp.specfit.modelerrs[1] for sp in self.speclist],unit=self.xarr.unit)
            self.fittable.add_column('widtherr',[sp.specfit.modelerrs[2] for sp in self.speclist],unit=self.xarr.unit)

    def ploteach(self, xunit=None, inherit_fit=False, plot_fit=True, plotfitkwargs={}, **plotkwargs):
        """
        Plot each spectrum in its own window
        inherit_fit - if specified, will grab the fitter & fitter properties from Spectra
        """
        for sp in self.speclist:
            if xunit is not None:
                sp.xarr.convert_to_unit(xunit,quiet=True)
            if inherit_fit:
                sp.specfit.fitter = self.specfit.fitter
                sp.specfit.modelpars = self.specfit.modelpars
                sp.specfit.model = np.interp(sp.xarr.as_unit(self.xarr.unit),
                                             self.xarr,self.specfit.fullmodel)

            sp.plotter(**plotkwargs)

            if plot_fit and self.specfit.model is not None:
                sp.specfit.plot_fit(**plotfitkwargs)


class ObsBlock(Spectra):
    """
    An Observation Block

    Consists of multiple spectra with a shared X-axis.  Intended to hold groups
    of observations of the same object in the same setup for later averaging.

    ObsBlocks can be indexed like python lists.
    """

    def __init__(self, speclist, xtype='frequency', xarr=None, force=False,
                 model_registry=None, **kwargs):

        if xarr is None:
            self.xarr = speclist[0].xarr
        else:
            self.xarr = xarr

        self.unit = speclist[0].unit
        self.header = speclist[0].header
        self.parse_header(self.header)

        for spec in speclist:
            if not isinstance(spec,Spectrum):
                raise TypeError("Must create an ObsBlock with a list of spectra.")
            if not np.array_equal(spec.xarr, self.xarr):
                if not force:
                    raise ValueError("Mismatch between X axes in ObsBlock")
            if spec.unit != self.unit:
                raise ValueError("Mismatched units")

        if force:
            self.speclist = [interpolation.interp(spec,self) for spec in speclist]
        else:
            self.speclist = speclist
        self.nobs = len(self.speclist)

        # Create a 2-dimensional array of the data
        self.data = np.array([sp.data for sp in self.speclist]).swapaxes(0,1).squeeze()
        self.error = np.array([sp.error for sp in self.speclist]).swapaxes(0,1).squeeze()

        self.plotter = plotters.Plotter(self)
        self._register_fitters(registry=model_registry)
        self.specfit = fitters.Specfit(self, Registry=self.Registry)
        self.baseline = baseline.Baseline(self)

    def average(self, weight=None, inverse_weight=False, error='erravgrtn', debug=False):
        """
        Average all scans in an ObsBlock.  Returns a single Spectrum object

        Parameters
        ----------
        weight : string
            a header keyword to weight by.   If not specified, the spectra will be
            averaged without weighting
        inverse_weight : bool
            Is the header keyword an inverse-weight (e.g., a variance?)
        error : ['scanrms','erravg','erravgrtn']
            estimate the error spectrum by one of three methods.
            'scanrms'   : the standard deviation of each pixel across all scans
            'erravg'    : the average of all input error spectra
            'erravgrtn' : the average of all input error spectra divided by sqrt(n_obs)
        """

        wtarr = np.isfinite(self.data)
        if weight is not None:
            if inverse_weight:
                for ii,sp in enumerate(self.speclist):
                    wtarr[:,ii] *= 1.0/sp.header.get(weight)
            else:
                for ii,sp in enumerate(self.speclist):
                    wtarr[:,ii] *= sp.header.get(weight)

        if self.header.get('EXPOSURE'):
            self.header['EXPOSURE'] = np.sum([sp.header['EXPOSURE'] for sp in self.speclist])

        data_nonan = np.nan_to_num(self.data)
        weighted_data = (data_nonan * wtarr)
        weighted_data_axsum = weighted_data.sum(axis=1)
        weight_axsum = wtarr.sum(axis=1)
        avgdata = weighted_data_axsum / weight_axsum
        if error == 'scanrms':
            # axis swapping is for projection... avgdata = 0'th axis
            errspec = np.sqrt((((data_nonan.swapaxes(0,1)-avgdata) *
                                wtarr.swapaxes(0,1))**2 /
                               wtarr.swapaxes(0,1)**2).swapaxes(0,1).sum(axis=1)
                             )
        elif error == 'erravg':
            errspec = self.error.mean(axis=1)
        elif error == 'erravgrtn':
            errspec = self.error.mean(axis=1) / np.sqrt(self.error.shape[1])

        spec = Spectrum(data=avgdata,error=errspec,xarr=self.xarr.copy(),header=self.header)
        spec._arithmetic_threshold = self._arithmetic_threshold

        if debug:
            # this statement, and much of the text above, is to test an absolutely insane error:
            # wtarr.sum(axis=1) randomly - say, one out of every 10-100 occurrences - fills with
            # nonsense values (1e-20, 1e-55, whatever).  There is no pattern to this; it occurs in
            # while loops, but ONLY IN THIS FUNCTION.  This is unreproduceable anywhere else.
            print("selfdata    min: %10g max: %10g" % (self.data.min(), self.data.max()))
            print("nonandata   min: %10g max: %10g" % (data_nonan.min(), data_nonan.max()))
            print("avgdata     min: %10g max: %10g" % (avgdata.min(), avgdata.max()))
            print("weight      sum: %10g" % (wtarr.sum()))
            print("data*weight sum: %10g" % ((data_nonan*wtarr).sum()))
            if np.abs(data_nonan.min()/avgdata.min()) > 1e10:
                import pdb; pdb.set_trace()

        return spec

    def __len__(self):
        return len(self.speclist)

    def __getitem__(self,index):
        """
        Can index Spectra to get the component Spectrum objects
        """
        return self.speclist[index]

    def smooth(self,smooth,**kwargs):
        """
        Smooth the spectrum by factor "smooth".  Options are defined in sm.smooth
        """
        smooth = round(smooth)
        self.data = sm.smooth_multispec(self.data,smooth,**kwargs)
        self.xarr = self.xarr[::smooth]
        if len(self.xarr) != len(self.data):
            raise ValueError("Convolution resulted in different X and Y array lengths.  Convmode should be 'same'.")
        self.error = sm.smooth_multispec(self.error,smooth,**kwargs)
        self.baseline.downsample(smooth)
        self.specfit.downsample(smooth)

        self._smooth_header(smooth)

class XCorrSpectrum(Spectrum):
    """ extraordinarily thin spectrum; just a name right now """
    pass

if __name__ == "__main__":
    import doctest
    doctest.testmod()
