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
import numpy as np
import smooth as sm
import pyfits
import readers,plotters,writers,baseline,units,measurements,speclines,arithmetic
import moments as moments_module
import fitters
import models
import history
import copy
from pyspeckit.specwarnings import warn
try:
    import atpy
    atpyOK = True
except ImportError:
    atpyOK = False

class Spectrum(object):
    """
    The core class for the spectroscopic toolkit.  Contains the data and error
    arrays along with wavelength / frequency / velocity information in various
    formats.

    Contains functions / classes to:
        -read and write FITS-compliant spectroscopic data sets
            -read fits binaries? (not implemented)
        -plot a spectrum
        -fit a spectrum both interactively and non-interactively
            -with gaussian, voigt, lorentzian, and multiple (gvl) profiles
        -fit a continuum "baseline" to a selected region with an n-order
          polynomial
        -perform fourier transforms and operations in fourier space on a spectrum (not implemented)

    Example usage:

    # automatically load a "compliant" (linear X-axis) FITS file

    spec = pyspeckit.Spectrum('test.fits')

    # plot the spectrum

    spec.plotter()


    """

    from arithmetic import interpnans

    def __init__(self, filename=None, filetype=None, xarr=None, data=None,
            error=None, header=None, doplot=False, maskdata=True,
            plotkwargs={}, xarrkwargs={}, **kwargs):
        """
        __init__
        Initialize the Spectrum.  Accepts files in the following formats:
            - .fits
            - .txt
            - .hdf5

        Must either pass in a filename or ALL of xarr, data, and header, plus
        optionally error.

        doplot - if specified, will generate a plot when created

        kwargs are passed to the reader, not the plotter
        """

        if filename:
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
            if filetype in ('fits','tspec','pyfits'):
                self.parse_header(self.header)
            elif filetype is 'txt':
                self.parse_text_header(self.header)
            elif filetype in ('hdf5', 'h5'):
                self.parse_hdf5_header(self.header)

            if isinstance(filename,str):
                self.fileprefix = filename.rsplit('.', 1)[0]    # Everything prior to .fits or .txt
        elif xarr is not None and data is not None:
            # technically, this is unpythonic.  But I don't want to search for all 10 attributes required.
            if issubclass(type(xarr),units.SpectroscopicAxis):
                self.xarr = xarr
            else:
                self.xarr = units.SpectroscopicAxis(xarr, **xarrkwargs)
            self.data = data
            if error is not None:
                self.error = error
            else:
                self.error = data * 0
            if hasattr(header,'get'):
                self.header = header
            else: # set as blank
                warn( "WARNING: Blank header." )
                self.header = pyfits.Header()
            self.parse_header(self.header)

        if maskdata:
            if hasattr(self.data,'mask'):
                self.data.mask += np.isnan(self.data) + np.isinf(self.data)
                self.error.mask += np.isnan(self.data) + np.isinf(self.data)
            else:
                self.data = np.ma.masked_where(np.isnan(self.data) + np.isinf(self.data), self.data)
                self.error = np.ma.masked_where(np.isnan(self.data) + np.isinf(self.data), self.error)

        self.plotter = plotters.Plotter(self)
        self._register_fitters()
        self.specfit = fitters.Specfit(self,Registry=self.Registry)
        self.baseline = baseline.Baseline(self)
        self.speclines = speclines
        self._sort()

        if doplot: self.plotter(**plotkwargs)

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

        for modelname, model in registry.multifitters.iteritems():
            self.Registry.add_fitter(modelname, model,
                    registry.npars[modelname], multisingle='multi',
                    key=registry.associated_keys.get(modelname))
        for modelname, model in registry.singlefitters.iteritems():
            self.Registry.add_fitter(modelname, model,
                    registry.npars[modelname], multisingle='single',
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
            self.xarr.dxarr = np.diff(self.xarr)
        
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
            self.xarr.units = Table.columns[xtype].unit
        elif xtype in units.xtype_dict:
            self.xarr.xtype = units.xtype_dict[xtype]
            self.xarr.units = Table.columns[xtype].unit
        else:
            print "Invalid xtype in text header - this may mean no text header was available.  X-axis units will be pixels"
            self.xarr.xtype = 'pixels'
            self.xarr.units = 'none'
            #raise ValueError("Invalid xtype in text header")
        self.ytype = Table.data.dtype.names[Table.datacol]
        self.units = Table.columns[self.ytype].unit
        self.header = pyfits.Header()
        self._update_header()

    def _update_header(self):
        self.header.update('CUNIT1',self.xarr.units)
        self.header.update('CTYPE1',self.xarr.xtype)
        self.header.update('BUNIT',self.units)
        self.header.update('BTYPE',self.ytype)
        
    def parse_hdf5_header(self, hdr):
        """
        HDF5 reader will create a hdr dictionary from HDF5 dataset
        attributes if they exist.  This routine will convert that dict
        to a pyfits header instance.
        """    
        
        self.xarr.xtype = hdr['xtype']
        self.xarr.xunits = hdr['xunits']
        self.ytype = hdr['ytype']
        self.units = hdr['yunits']
        self.header = pyfits.Header()
        self.header.update('CUNIT1', self.xarr.xunits)
        self.header.update('CTYPE1', self.xarr.xtype)
        self.header.update('BUNIT', self.ytype)
        self.header.update('BTYPE', self.units)

    def parse_header(self,hdr,specname=None):
        """
        Parse parameters from a .fits header into required spectrum structure
        parameters

        This should be moved to the FITSSpectrum subclass when that is available
        """

        if hdr.get('BUNIT'):
            self.units = hdr.get('BUNIT').strip()
        else:
            self.units = 'undefined'
            
        if hdr.get('BTYPE'):
            self.ytype = hdr.get('BTYPE').strip()
        else:
            self.ytype = 'data'

        if specname is not None:
            self.specname = specname
        elif hdr.get('OBJECT'):
            self.specname = hdr.get('OBJECT')
        else:
            self.specname = ''
            
    def measure(self, z=None, d=None, fluxnorm=None, miscline=None, misctol=None, ignore=None, derive=True, **kwargs):
        """
        Initialize the measurements class - only do this after you have run a fitter otherwise pyspeckit will be angry!
        """
        self.measurements=measurements.Measurements(self, z=z, d=d,
                fluxnorm=fluxnorm, miscline=miscline, misctol=misctol,
                ignore=ignore, derive=derive, **kwargs)

    def crop(self, x1, x2, units=None, **kwargs):
        """
        Replace the current spectrum with a subset from x1 to x2 in current
        units

        Fixes CRPIX1 and baseline and model spectra to match cropped data spectrum

        """
        # do slice (this code is redundant... need to figure out how to fix that)
        x1pix = np.argmin(np.abs(x1-self.xarr.as_unit(units)))
        x2pix = np.argmin(np.abs(x2-self.xarr.as_unit(units)))
        if x1pix > x2pix: x1pix,x2pix = x2pix,x1pix
        if x1pix == x2pix:
            raise IndexError("ERROR: Trying to crop to zero size.")

        self = self.slice(x1, x2, copy=False, **kwargs)
        # a baseline spectrum is always defined, even if it is all zeros
        # this is needed to prevent size mismatches.  There may be a more
        # elegant way to do this...
        self.baseline.crop(x1pix,x2pix)
        self.specfit.crop(x1pix,x2pix)

        if hasattr(self,'header'):
            history.write_history(self.header,"CROP: Cropped from %g to %g (pixel %i to %i)" % (x1,x2,x1pix,x2pix))

            if self.header.get('CRPIX1'):
                self.header.update('CRPIX1',self.header.get('CRPIX1') - x1pix)
                history.write_history(self.header,"CROP: Changed CRPIX1 from %f to %f" % (self.header.get('CRPIX1')+x1pix,self.header.get('CRPIX1')))

    def slice(self, x1, x2, units=None, copy=True):
        """
        Slice a spectrum.  Defaults to the current units.  Can specify 'pixels'
        if you want to just slice by index
        """
        # do slice
        x1pix = np.argmin(np.abs(x1-self.xarr.as_unit(units)))
        x2pix = np.argmin(np.abs(x2-self.xarr.as_unit(units)))
        if x1pix > x2pix: x1pix,x2pix = x2pix,x1pix
        if x1pix == x2pix:
            raise IndexError("ERROR: Trying to crop to zero size.")

        if copy:
            sp = self.copy()
        else:
            sp = self
        # reset the plot window when cropping
        sp.plotter.xmin = sp.xarr[x1pix]
        sp.plotter.xmax = sp.xarr[x2pix]

        sp.xarr = sp.xarr[x1pix:x2pix]
        sp.data = sp.data[x1pix:x2pix]
        sp.error = sp.error[x1pix:x2pix]

        return sp


    def smooth(self,smooth,**kwargs):
        """
        Smooth the spectrum by factor "smooth".  
        """
        smooth = round(smooth)
        self.data = sm.smooth(self.data,smooth,**kwargs)
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
            self.header.update('CDELT1',self.header.get('CDELT1') * float(smooth))
            self.header.update('CRPIX1',self.header.get('CRPIX1') / float(smooth))

            history.write_history(self.header,"SMOOTH: Smoothed and downsampled spectrum by factor %i" % (smooth))
            history.write_history(self.header,"SMOOTH: Changed CRPIX1 from %f to %f" % (self.header.get('CRPIX1')*float(smooth),self.header.get('CRPIX1')))
            history.write_history(self.header,"SMOOTH: Changed CDELT1 from %f to %f" % (self.header.get('CRPIX1')/float(smooth),self.header.get('CRPIX1')))

    def shape(self):
        """
        Return the data shape
        """
        return self.data.shape

    def __len__(self):
        return len(self.data)

    def __repr__(self):
        return r'<Spectrum object over spectral range %6.5g : %6.5g %s and flux range = [%2.1f, %2.1f] %s>' % \
                (self.xarr.min(), self.xarr.max(), self.xarr.units, self.data.min(), self.data.max(), self.units)
    

    def copy(self):
        """
        Create a copy of the spectrum with its own plotter, fitter, etc.
        Useful for, e.g., comparing smoothed to unsmoothed data
        """

        newspec = copy.copy(self)
        newspec.plotter = plotters.Plotter(newspec)
        newspec._register_fitters()
        newspec.specfit = fitters.Specfit(newspec,Registry=newspec.Registry)
        newspec.baseline = baseline.Baseline(newspec)

        return newspec

    def stats(self, statrange=(), interactive=False):
        """
        Return some statistical measures in a dictionary (somewhat self-explanatory)

        range - X-range over which to perform measures
        interactive - specify range interactively in plotter

        """

        if len(statrange) == 2:
            pix1 = self.xarr.x_to_pix(statrange[0])
            pix2 = self.xarr.x_to_pix(statrange[1])
            if pix1 > pix2: pix1,pix2 = pix2,pix1
            elif pix1 == pix2: raise ValueError("Invalid statistics range - includes 0 pixels")
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
        if not linetype in self.__dict__: # don't replace it if it already exists
            self.__dict__[linetype] = speclines.__dict__[linetype].__dict__[linetype+"_lines"](self,**kwargs)

    def moments(self, unit='km/s', **kwargs):
        """
        Return the moments of the spectrum.  In order to assure that the 1st
        and 2nd moments are meaningful, a 'default' unit is set.  If unit is not
        set, will use current unit.

        """

        if unit is False or unit is None:
            return moments_module.moments(self.xarr, self.data, **kwargs)
        else:
            return moments_module.moments(self.xarr.as_unit(unit), self.data, **kwargs)

    moments.__doc__ += moments_module.moments.__doc__



class Spectra(Spectrum):
    """
    A list of individual Spectrum objects.  Intended to be used for
    concatenating different wavelength observations of the SAME OBJECT.  Can be
    operated on just like any Spectrum object, incuding fitting.  Useful for
    fitting multiple lines on non-continguous axes simultaneously.  Be wary of
    plotting these though...

    Can be indexed like python lists.  

    X array is forcibly sorted in increasing order
    """

    def __init__(self, speclist, xunits='GHz', **kwargs):
        print "Creating spectra"
        speclist = list(speclist)
        for ii,spec in enumerate(speclist):
            if type(spec) is str:
                spec = Spectrum(spec)
                speclist[ii] = spec

        self.speclist = speclist

        print "Concatenating data"
        self.xarr = units.SpectroscopicAxes([sp.xarr.as_unit(xunits) for sp in speclist])
        self.xarr.units = xunits 
        self.xarr.xtype = units.unit_type_dict[xunits]
        self.data = np.ma.concatenate([sp.data for sp in speclist])
        self.error = np.ma.concatenate([sp.error for sp in speclist])
        self._sort()

        self.header = pyfits.Header()
        for spec in speclist:
            for key,value in spec.header.items():
                self.header.update(key,value)

        self.plotter = plotters.Plotter(self)
        self._register_fitters()
        self.specfit = fitters.Specfit(self,Registry=self.Registry)
        self.baseline = baseline.Baseline(self)
        
        self.units = speclist[0].units
        for spec in speclist: 
            if spec.units != self.units: 
                raise ValueError("Mismatched units")

    def __add__(self,other):
        """
        Adding "Spectra" together will concatenate them
        """
        if type(other) is Spectra or Spectrum:
            self.speclist += other.speclist
        elif type(other) is Spectrum:
            self.speclist += [other.speclist]
        if other.units != self.units:
            raise ValueError("Mismatched units")

        if other.xarr.units != self.xarr.units:
            # convert all inputs to same units
            spec.xarr.convert_to_units(self.xarr.units,**kwargs)
        self.xarr = units.SpectroscopicAxes([self.xarr,spec.xarr])
        self.data = np.concatenate([self.data,spec.data])
        self.error = np.concatenate([self.error,spec.error])
        self._sort()

    def __getitem__(self,index):
        """
        Can index Spectra to get the component Spectrum objects
        """
        return self.speclist[index]

    def __len__(self): return len(self.speclist)

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
            self.fittable.add_column('amplitude',[sp.specfit.modelpars[0] for sp in self.speclist],unit=self.units)
            self.fittable.add_column('center',[sp.specfit.modelpars[1] for sp in self.speclist],unit=self.xarr.units)
            self.fittable.add_column('width',[sp.specfit.modelpars[2] for sp in self.speclist],unit=self.xarr.units)
            self.fittable.add_column('amplitudeerr',[sp.specfit.modelerrs[0] for sp in self.speclist],unit=self.units)
            self.fittable.add_column('centererr',[sp.specfit.modelerrs[1] for sp in self.speclist],unit=self.xarr.units)
            self.fittable.add_column('widtherr',[sp.specfit.modelerrs[2] for sp in self.speclist],unit=self.xarr.units)

    def ploteach(self, xunits=None, inherit_fit=False, plot_fit=True, plotfitkwargs={}, **plotkwargs):
        """
        Plot each spectrum in its own window
        inherit_fit - if specified, will grab the fitter & fitter properties from Spectra
        """
        for sp in self.speclist:
            if xunits is not None:
                sp.xarr.convert_to_unit(xunits,quiet=True)
            if inherit_fit:
                sp.specfit.fitter = self.specfit.fitter
                sp.specfit.modelpars = self.specfit.modelpars
                sp.specfit.model = np.interp(sp.xarr.as_unit(self.xarr.units),self.xarr,self.specfit.fullmodel)

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

    def __init__(self, speclist, xtype='frequency', xarr=None, force=False, **kwargs):

        if xarr is None:
            self.xarr = speclist[0].xarr
        else:
            self.xarr = xarr

        self.units = speclist[0].units
        self.header = speclist[0].header
        self.parse_header(self.header)

        for spec in speclist:
            if type(spec) is not Spectrum:
                raise TypeError("Must create an ObsBlock with a list of spectra.")
            if not np.array_equal(spec.xarr, self.xarr):
                if not force:
                    raise ValueError("Mismatch between X axes in ObsBlock")
            if spec.units != self.units: 
                raise ValueError("Mismatched units")

        if force:
            self.speclist = [arithmetic.interp(spec,self) for spec in speclist]
        else:
            self.speclist = speclist
        self.nobs = len(self.speclist)

        # Create a 2-dimensional array of the data
        self.data = np.array([sp.data for sp in self.speclist]).swapaxes(0,1)
        self.error = np.array([sp.error for sp in self.speclist]).swapaxes(0,1)

        self.plotter = plotters.Plotter(self)
        self._register_fitters()
        self.specfit = fitters.Specfit(self,Registry=self.Registry)
        self.baseline = baseline.Baseline(self)
        
    def average(self, weight=None, inverse_weight=False, error='erravgrtn'):
        """
        Average all scans in an ObsBlock.  Returns a single Spectrum object

        weight - a header keyword to weight by
        error - estimate the error spectrum.  Can be:
            'scanrms'   - the standard deviation of each pixel across all scans
            'erravg'    - the average of all input error spectra
            'erravgrtn' - the average of all input error spectra divided by sqrt(n_obs)
        """

        if weight is not None:
            if inverse_weight:
                wtarr = np.reshape( np.array([1.0/sp.header.get(weight) for sp in self.speclist]) , [self.nobs,1] )
            else:
                wtarr = np.reshape( np.array([sp.header.get(weight) for sp in self.speclist]) , [self.nobs,1] )
        else:
            wtarr = np.ones(self.data.shape[1])

        if self.header.get('EXPOSURE'):
            self.header['EXPOSURE'] = np.sum([sp.header['EXPOSURE'] for sp in self.speclist])

        avgdata = (self.data * wtarr.swapaxes(0,1)).sum(axis=1) / wtarr.sum()
        if error is 'scanrms':
            errspec = np.sqrt( (((self.data.swapaxes(0,1)-avgdata) * wtarr)**2 / wtarr**2).swapaxes(0,1).sum(axis=1) )
        elif error is 'erravg':
            errspec = self.error.mean(axis=1)
        elif error is 'erravgrtn':
            errspec = self.error.mean(axis=1) / np.sqrt(self.error.shape[1])

        spec = Spectrum(data=avgdata,error=errspec,xarr=self.xarr.copy(),header=self.header)

        return spec

    def __len__(self): return len(self.speclist)

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
