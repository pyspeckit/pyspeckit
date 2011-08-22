"""
The spectrum module consists of the Spectrum class, with child classes ObsBlock
and Spectra for multi-spectrum analysis of different types.

The Spectrum class is the main functional code.
ObsBlocks are containers of multiple spectra of different objects
The Spectra class is a container of multiple spectra of the *same* object at
   different wavelengths/frequencies

.. moduleauthor:: Adam Ginsburg <adam.g.ginsburg@gmail.com>
"""
import numpy as np
import smooth as sm
import pyfits
import readers,plotters,writers,baseline,units,measurements,speclines,arithmetic
import fitters
import models
import history
import copy
try:
    import atpy
    atpyOK = True
except ImportError:
    atpyOK = False

def register_fitter(Registry, name, function, npars, multisingle='single',
        override=False, key=None):
    ''' 
    Register a fitter function.

    Required Arguments:

        *name*: [ string ]
            The fit function name. 

        *function*: [ function ]
            The fitter function.  Single-fitters should take npars + 1 input
            parameters, where the +1 is for a 0th order baseline fit.  They
            should accept an X-axis and data and standard fitting-function
            inputs (see, e.g., gaussfitter).  Multi-fitters should take N *
            npars, but should also operate on X-axis and data arguments.

        *npars*: [ int ]
            How many parameters does the function being fit accept?

    Optional Keyword Arguments:

        *multisingle*: [ 'multi' | 'single' ] 
            Is the function a single-function fitter (with a background), or
            does it allow N copies of the fitting function?

        *override*: [ True | False ]
            Whether to override any existing type if already present.

        *key*: [ char ]
            Key to select the fitter in interactive mode
    '''

    if multisingle == 'single':
        if not name in Registry.singlefitters or override:
            Registry.singlefitters[name] = function
    elif multisingle == 'multi':
        if not name in Registry.multifitters or override:
            Registry.multifitters[name] = function
    elif name in Registry.singlefitters or name in Registry.multifitters:
        raise Exception("Fitting function %s is already defined" % name)

    if key is not None:
        Registry.fitkeys[key] = name
        Registry.interactive_help_message += "\n'%s' - select fitter %s" % (key,name)
    Registry.npars[name] = npars



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

    """

    def __init__(self,filename=None, filetype=None, xarr=None, data=None,
            error=None, header=None, doplot=False, plotkwargs={}, **kwargs):
        """
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
                if readers.suffix_types.has_key(suffix):
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

            if isinstance(filename,str):
                self.fileprefix = filename.rsplit('.', 1)[0]    # Everything prior to .fits or .txt
        elif xarr is not None and data is not None:
            self.xarr = xarr
            self.data = data
            if error is not None:
                self.error = error
            else:
                self.error = data * 0
            self.header = header
            self.parse_header(header)


        self.plotter = plotters.Plotter(self)
        self._register_fitters()
        self.specfit = fitters.Specfit(self,Registry=self.Registry)
        self.baseline = baseline.Baseline(self)
        self.speclines = speclines

        if doplot: self.plotter(**plotkwargs)

    def _register_fitters(self):
        """
        Register fitters independently for each spectrum instance
        """
        Registry = fitters.Registry()
        register_fitter(Registry,'ammonia',models.ammonia_model(multisingle='multi'),6,multisingle='multi',key='a')
        # not implemented register_fitter(Registry,'ammonia',models.ammonia_model(multisingle='single'),6,multisingle='single',key='A')
        register_fitter(Registry,'formaldehyde',models.formaldehyde_model(multisingle='multi'),3,multisingle='multi',key='F') # CAN'T USE f!  reserved for fitting
        register_fitter(Registry,'formaldehyde',models.formaldehyde_model(multisingle='single'),3,multisingle='single')
        register_fitter(Registry,'gaussian',models.gaussian_fitter(multisingle='multi'),3,multisingle='multi',key='g')
        register_fitter(Registry,'gaussian',models.gaussian_fitter(multisingle='single'),3,multisingle='single')
        register_fitter(Registry,'voigt',models.voigt_fitter(multisingle='multi'),4,multisingle='multi',key='v')
        register_fitter(Registry,'voigt',models.voigt_fitter(multisingle='single'),4,multisingle='single')
        self.Registry = Registry

        
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
            if writers.suffix_types.has_key(suffix):
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
        self.xarr.xtype = Table.data.dtype.names[0]
        self.xarr.xunits = Table.columns[self.xarr.xtype].unit
        self.ytype = Table.data.dtype.names[1]
        self.units = Table.columns[self.ytype].unit
        self.header = pyfits.Header()
        self.header.update('CUNIT1',self.xarr.xunits)
        self.header.update('CTYPE1',self.xarr.xtype)
        self.header.update('BUNIT',self.units)
        self.header.update('BTYPE',self.ytype)

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

        if specname is not None:
            self.specname = specname
        elif hdr.get('OBJECT'):
            self.specname = hdr.get('OBJECT')
        else:
            self.specname = ''
            
    def measure(self, z = None, d = None, fluxnorm = None, miscline = None, misctol = None, ignore = None):
        """
        Initialize the measurements class - only do this after you have run a fitter otherwise pyspeckit will be angry!
        """
        self.measurements = measurements.Measurements(self, z = z, d = d, fluxnorm = fluxnorm, miscline = miscline, misctol = misctol, ignore = ignore)                

    def crop(self,x1,x2):
        """
        Replace the current spectrum with a subset from x1 to x2 in current
        units

        Fixes CRPIX1 and baseline and model spectra to match cropped data spectrum

        """
        x1pix = np.argmin(np.abs(x1-self.xarr))
        x2pix = np.argmin(np.abs(x2-self.xarr))
        if x1pix > x2pix: x1pix,x2pix = x2pix,x1pix

        self.plotter.xmin = self.xarr[x1pix]
        self.plotter.xmax = self.xarr[x2pix]
        self.xarr = self.xarr[x1pix:x2pix]
        self.data = self.data[x1pix:x2pix]
        self.error = self.error[x1pix:x2pix]
        # a baseline spectrum is always defined, even if it is all zeros
        # this is needed to prevent size mismatches.  There may be a more
        # elegant way to do this...
        self.baseline.crop(x1pix,x2pix)
        self.specfit.crop(x1pix,x2pix)

        if self.header.get('CRPIX1'):
            self.header.update('CRPIX1',self.header.get('CRPIX1') - x1pix)

        history.write_history(self.header,"CROP: Cropped from %g to %g (pixel %i to %i)" % (x1,x2,x1pix,x2pix))
        history.write_history(self.header,"CROP: Changed CRPIX1 from %f to %f" % (self.header.get('CRPIX1')+x1pix,self.header.get('CRPIX1')))

    def smooth(self,smooth,**kwargs):
        """
        Smooth the spectrum by factor "smooth".  Options are defined in sm.smooth
        """
        smooth = round(smooth)
        self.data = sm.smooth(self.data,smooth,**kwargs)
        self.xarr = self.xarr[::smooth]
        if len(self.xarr) != len(self.data):
            raise ValueError("Convolution resulted in different X and Y array lengths.  Convmode should be 'same'.")
        self.error = sm.smooth(self.error,smooth,**kwargs)
        self.baseline.downsample(smooth)
        self.specfit.downsample(smooth)
    
        self._smooth_header(smooth)

    def _smooth_header(self,smooth):
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

    def stats(self):
        """
        Return some statistical measures
        """

        stats = {
            "npts": self.data.shape[0],
            "std": self.data.std(),
            "mean": self.data.mean(),
            "median": np.median(self.data),
            "min": self.data.min(),
            "max": self.data.max(),}
        return stats

    def getlines(self, linetype='radio'):
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
        if not self.__dict__.has_key(linetype): # don't replace it if it already exists
            self.__dict__[linetype] = speclines.__dict__[linetype].__dict__[linetype+"_lines"](self)




class Spectra(Spectrum):
    """
    A list of individual Spectrum objects.  Intended to be used for
    concatenating different wavelength observations of the SAME OBJECT.  Can be
    operated on just like any Spectrum object, incuding fitting.  Useful for
    fitting multiple lines on non-continguous axes simultaneously.  Be wary of
    plotting these though...
    """

    def __init__(self,speclist,xunits='GHz',**kwargs):
        print "Creating spectra"
        speclist = list(speclist)
        for ii,spec in enumerate(speclist):
            if type(spec) is str:
                spec = Spectrum(spec)
                speclist[ii] = spec
            if spec.xarr.units != xunits:
                # convert all inputs to same units
                spec.xarr.convert_to_unit(xunits, **kwargs)

        self.speclist = speclist

        print "Concatenating data"
        self.xarr = units.SpectroscopicAxes([sp.xarr for sp in speclist])
        self.data = np.concatenate([sp.data for sp in speclist])
        self.error = np.concatenate([sp.error for sp in speclist])
        self._sort()

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

    def _sort(self):
        """ Sort the data in order of increasing X (could be decreasing, but
        must be monotonic for plotting reasons) """

        indices = self.xarr.argsort()
        self.xarr = self.xarr[indices]
        self.data = self.data[indices]
        self.error = self.error[indices]

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


class ObsBlock(Spectra):
    """
    An Observation Block

    Consists of multiple spectra with a shared X-axis.  Intended to hold groups
    of observations of the same object in the same setup for later averaging.
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
            if not (spec.xarr == self.xarr).all():
                if force:
                    spec = arithmetic.interp(spec,self)
                else:
                    raise ValueError("Mismatch between X axes in ObsBlock")
            if spec.units != self.units: 
                raise ValueError("Mismatched units")

        self.speclist = speclist
        self.nobs = len(speclist)

        # Create a 2-dimensional array of the data
        self.data = np.array([sp.data for sp in speclist]).swapaxes(0,1)
        self.error = np.array([sp.error for sp in speclist]).swapaxes(0,1)

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
