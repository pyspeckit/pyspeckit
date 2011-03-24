import numpy as np
import smooth as sm
import pyfits
import readers,fitters,plotters,writers,baseline,units,speclines


class Spectrum(object):
    """
    The core class for the spectroscopic toolkit.  Contains the data and error
    arrays along with wavelength / frequency / velocity information in various
    formats.

    Will contain functions / classes to:
        -read and write FITS-compliant spectroscopic data sets
            -read fits binaries?
        -plot a spectrum
        -fit a spectrum both interactively and non-interactively
            -with gaussian, voigt, lorentzian, and multiple (gvl) profiles
        -fit a continuum "baseline" to a selected region with an n-order
          polynomial
        -perform fourier transforms and operations in fourier space on a spectrum
    """

    def __init__(self,filename, filetype=None, **kwargs):
        """
        Initialize the Spectrum.  Accepts files in the following formats:
            - .fits
            - .txt (not functional yet)

        Ideally, we'd like to do something where each independent type of data
        file has its own subclass for file reading, but we haven't figured out
        how to automatically determine which subclass to use yet.
        """

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
        self.plotter = plotters.Plotter(self)
        self.specfit = fitters.Specfit(self)
        self.baseline = baseline.Baseline(self)
        
        # Initialize writers
        self.writer = {}
        for writer in writers.writers: self.writer[writer] = writers.writers[writer](self)

    def parse_text_header(self,Table):
        """
        Grab relevant parameters from a table header (xaxis type, etc)
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

    def crop(self,x1,x2):
        """
        Replace the current spectrum with a subset from x1 to x2 in current
        units
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

    def smooth(self,smooth,**kwargs):
        """
        Smooth the spectrum.  Options are defined in sm.smooth
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


class Spectra(Spectrum):
    """
    A list of individual Spectrum objects.  Can be operated on just like any
    Spectrum object, incuding fitting.  Useful for fitting multiple lines on
    non-continguous axes simultaneously.  Be wary of plotting these though...
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
        self.writer = writers.Writer(self)

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
        return self.speclist[index]
