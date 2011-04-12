import numpy as np
import smooth as sm
import pyfits
import readers,fitters,plotters,writers,baseline,units,speclines
import history


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
            error=None, header=None, **kwargs):
        """
        Initialize the Spectrum.  Accepts files in the following formats:
            - .fits
            - .txt
            - hdf5

        Must either pass in a filename or ALL of xarr, data, and header, plus
        optionally error.

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

            self.data,self.error,self.xarr,self.header = reader(filename)    
            
            # these should probably be replaced with registerable function s...
            if filetype in ('fits','tspec'):
                self.parse_header(self.header)
            elif filetype is 'txt':
                self.parse_text_header(self.header)

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
        self.specfit = fitters.Specfit(self)
        self.baseline = baseline.Baseline(self)
        
        # Initialize writers
        self.writer = {}
        for writer in writers.writers: self.writer[writer] = writers.writers[writer](self)

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

    def parse_header(self,hdr,specname=None, wcstype=''):
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


class Spectra(Spectrum):
    """
    A list of individual Spectrum objects.  Intended to be used for
    concatenating different wavelength observations of the SAME OBJECT.  Can be
    operated on just like any Spectrum object, incuding fitting.  Useful for
    fitting multiple lines on non-continguous axes simultaneously.  Be wary of
    plotting these though...
    """

    def __init__(self,speclist,xtype='frequency',**kwargs):
        print "Creating spectra"
        for spec in speclist:
            if type(spec) is str:
                spec = Spectrum(spec)
            if spec.xarr.xtype is not xtype:
                # convert all inputs to same units
                spec.xarr.change_xtype(xtype,**kwargs)

        self.speclist = speclist

        print "Concatenating data"
        self.xarr = units.SpectroscopicAxes([sp.xarr for sp in speclist])
        self.data = np.concatenate([sp.data for sp in speclist])
        self.error = np.concatenate([sp.error for sp in speclist])

        self.plotter = plotters.Plotter(self)
        self.specfit = fitters.Specfit(self)
        self.baseline = baseline.Baseline(self)
        
        # Initialize writers
        self.writer = {}
        for writer in writers.writers: self.writer[writer] = writers.writers[writer](self)

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

        if other.xarr.xtype is not self.xarr.xtype:
            # convert all inputs to same units
            spec.xarr.change_xtype(self.xarr.xtype,**kwargs)
        self.xarr = units.SpectroscopicAxes([self.xarr,spec.xarr])
        self.data = np.concatenate([self.data,spec.data])
        self.error = np.concatenate([self.error,spec.error])

    def __getitem__(self,index):
        """
        Can index Spectra to get the component Spectrum objects
        """
        return self.speclist[index]

    def __len__(self): return len(self.speclist)

class ObsBlock(Spectrum):
    """
    An Observation Block

    Consists of multiple spectra with a shared X-axis.  Intended to hold groups
    of observations of the same object in the same setup for later averaging.
    """

    def __init__(self,speclist,xtype='frequency',xarr=None,**kwargs):

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
                raise ValueError("Mismatch between X axes in ObsBlock")
            if spec.units != self.units: 
                raise ValueError("Mismatched units")

        self.speclist = speclist

        # Create a 2-dimensional array of the data
        self.data = np.array([sp.data for sp in speclist]).swapaxes(0,1)
        self.error = np.array([sp.error for sp in speclist]).swapaxes(0,1)

        self.plotter = plotters.Plotter(self)
        self.specfit = fitters.Specfit(self)
        self.baseline = baseline.Baseline(self)
        
        # Initialize writers
        self.writer = {}
        for writer in writers.writers: self.writer[writer] = writers.writers[writer](self)

    def average(self, weight=None, error='erravgrtn'):
        """
        Average all scans in an ObsBlock.  Returns a single Spectrum object

        weight - a header keyword to weight by
        error - estimate the error spectrum.  Can be:
            'scanrms'   - the standard deviation of each pixel across all scans
            'erravg'    - the average of all input error spectra
            'erravgrtn' - the average of all input error spectra divided by sqrt(n_obs)
        """

        if weight is not None:
            wtarr = np.array([sp.header.get(weight) for sp in self.speclist])
        else:
            wtarr = np.ones(self.data.shape[1])

        if self.header.get('EXPOSURE'):
            self.header['EXPOSURE'] = np.sum([sp.header['EXPOSURE'] for sp in self.speclist])

        avgdata = (self.data * wtarr).sum(axis=1) / wtarr.sum()
        if error is 'scanrms':
            errspec = np.sqrt( (((self.data-avgdata) * wtarr)**2 / wtarr**2).sum(axis=1) )
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
